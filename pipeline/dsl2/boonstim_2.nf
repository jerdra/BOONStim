nextflow.preview.dsl=2

if (!params.bids || !params.out){

    log.info('Insufficient input specification')
    log.info('Needs --bids, --out!')
    log.info('Exiting')
    System.exit(1)

}

log.info("BIDS Directory: $params.bids")
log.info("Output Directory: $params.out")

if (!params.anat_invocation || !params.ciftify_invocation || !params.ciftify_descriptor \
    || !params.anat_descriptor){
    log.info("Missing BOSH invocations and descriptor for fmriprep and ciftify!")
    log.info("Exiting with Error")
    System.exit(1)
}

if (!params.subject){
    log.info("Missing subject argument!")
    log.info("Exiting...")
    System.exit(1)
}

log.info("Using Descriptor Files: $params.anat_descriptor and $params.ciftify_descriptor")
log.info("Using Invocation Files: $params.anat_invocation and $params.ciftify_invocation")
log.info("Using containers: $params.fmriprep_img and $params.ciftify_img")
log.info("Subject to run: $params.subject")

include cifti_meshing from './modules/cifti_mesh_2_wf.nf' params(params)
include weightfunc_wf from "${params.weightworkflow}" params(params)
include centroid_wf from './modules/centroid_wf.nf' params(params)
include parameterization_wf from './modules/surfparams_wf.nf' params(params)
include tet_project_wf from './modules/tetrahedral_wf.nf' params(params)

process optimize_coil {

    label 'rtms'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh), path(weights), path(C), path(R), path(bounds), path(coil)

    output:
    tuple val(sub), path('coil_position'), emit: position
    tuple val(sub), path('coil_orientation'), emit: orientation
    tuple val(sub), path('history'), emit: history

    shell:
    '''
    /scripts/optimize_fem.py !{msh} !{weights} !{C} !{bounds} !{R} !{coil} \
                             coil_position coil_orientation \
                             --history history \
                             --cpus !{params.cpus}
    '''
}

def surf_branch = branchCriteria {
                white: it[2].contains('white.32k_fs_LR')
                    return it
                pial: it[2].contains('pial.32k_fs_LR')
                    return it
                midthick: it[2].contains('midthickness.32k_fs_LR')
                    return it
            }

process construct_outputs{

    input:
        tuple val(sub), \
        path(msh), path(t1fs), path(m2m), \
        path(wfunc), path(mask), \
        path(loc), path(rot), path(hist)

    output:
        tuple val(sub), path("$sub"), emit: subject

    shell:
    '''
    #!/bin/bash

    mkdir optimization
    mv !{loc} !{rot} !{hist} optimization

    mkdir !{sub}
    mv * !{sub} || true
    '''


}

workflow {
    main:

        //Main preprocessing
        cifti_meshing(params.subject)

        //Weightfunc calculation
        weightfunc_input = cifti_meshing.fmriprep.join(cifti_meshing.cifti)
        weightfunc_wf(weightfunc_input)

        //Parameterization workflow
        cifti_meshing.cifti
                     .spread(['pial','white','midthickness'])
                     .spread(['L','R'])
                     .map{n,c,s,h -> [n,h,"${c}/T1w/fsaverage_LR32k/${n}.${h}.${s}.32k_fs_LR.surf.gii"]}
                     .branch(surf_branch).set { surfs }

        centroid_wf(weightfunc_wf.out.mask, \
                    surfs.pial, surfs.white, surfs.midthick, \
                    cifti_meshing.t1fs_conform)

        parameterization_wf(cifti_meshing.msh, \
                            centroid_wf.out.centroid)

        //FEM projection workflow
        tet_project_wf(weightfunc_wf.weightfunc, \
                       surfs.pial, surfs.white, surfs.midthickness, \
                       cifti_meshing.t1fs_conform, \
                       cifti_meshing.msh)

        //Optimization
        optimize_inputs = cifti_meshing.msh
                                       .join(tet_project_wf.out.fem_weights)
                                       .join(parameterization.out.C)
                                       .join(parameterization.out.R)
                                       .join(parameterization.out.bounds)
                                       .map{s,m,w,C,R,b -> [s,m,w,C,R,b,"$params.coil"]}
        optimize_coil(optimize_inputs)

        //Output
        construct_outputs_input = cifti_meshing.msh
                                               .join(cifti_meshing.t1fs_conform)
                                               .join(cifti_meshing.mesh_m2m)
                                               .join(cifti_meshing.mesh_fs)
                                               .join(weightfunc_wf.out.weightfunc)
                                               .join(weightfunc_wf.out.mask)
                                               .join(optimize_coil.out.position)
                                               .join(optimize_coil.out.orientation)
                                               .join(optimize_coil.out.history)
        construct_outputs(construct_outputs_input)

    publish:
        cifti_meshing.cifti to: "$params.out/ciftify", mode: 'copy'
        cifti_meshing.fmriprep to: "$params.out/fmriprep", mode: 'copy'
        construct_outputs.out.subject to: "$params.out/boonstim", mode: 'copy'
}

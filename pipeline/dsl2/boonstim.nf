// BOONSTIM FULL PIPELINE WORKFLOW
nextflow.preview.dsl=2

if (!params.bids || !params.out){

    log.info('Insufficient input specification!')
    log.info('Needs --bids, --out!')
    log.info('Exiting...')
    System.exit(1)

}

log.info("BIDS Directory: $params.bids")
log.info("Output Directory: $params.out")

//Now check BOSH invocation of Ciftify
if (!params.anat_invocation || !params.ciftify_invocation || !params.ciftify_descriptor || !params.anat_descriptor) {

    log.info('Missing BOSH invocations and descriptor for fmriprep-ciftify!')
    log.info('Please have an --anat-only invocation for fmriprep and an associated descriptor')
    log.info('Exiting with Error')
    System.exit(1)

}

//Subject list flag
if (params.subjects) {
    log.info ("Subject list file provided: $params.subjects")
}

log.info("Using Descriptor Files: $params.anat_descriptor and $params.ciftify_descriptor")
log.info("Using Invocation Files: $params.anat_invocation and $params.ciftify_invocation")
log.info("Using containers: $params.fmriprep_img and $params.ciftify_img")
log.info("Using user-defined ROI workflow: $params.weightworkflow")
log.info("Using Invocation Files: $params.anat_invocation and $params.ciftify_invocation")

///////////////////////////////////////////////////////////////////////////////

// IMPORT MODULES WORKFLOWS
include cifti_meshing from './modules/cifti_mesh_wf.nf' params(params)
include make_giftis from './modules/fs2gifti.nf' params(params)
include registration_wf from './modules/register_fs2cifti_wf.nf' params(params)
include weightfunc_wf from "${params.weightworkflow}" params(params)
include resample2native_wf as resamplemask_wf from './modules/resample2native.nf' params(params)
include resample2native_wf as resampleweightfunc_wf from './modules/resample2native.nf' params(params)
include centroid_wf from './modules/centroid_wf.nf' params(params)
include parameterization_wf from './modules/surfparams_wf.nf' params(params)
include tet_project_wf from './modules/tetrahedral_wf.nf' params(params)

// IMPORT MODULES PROCESSES
include apply_mask as centroid_mask from './modules/utils.nf' params(params)
include apply_mask as weightfunc_mask from './modules/utils.nf' params(params)
include cifti_dilate as dilate_mask from './modules/utils.nf' params(params)

///////////////////////////////////////////////////////////////////////////////

//// Extract subjects to run
all_dirs = file(params.bids).list()
input_dirs = new File(params.bids).list()
output_dirs = new File(params.out).list()

if (params.subjects) {
    sublist = file(params.subjects)
    bids_channel = Channel.from(sublist)
                               .splitText() { it.strip() }
} else{
    bids_channel = Channel.from(all_dirs)
                          .filter { it.contains('sub-') }
}

// Process definitions
process optimize_coil {

    stageInMode 'copy'
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
                             --n-iters 30 \
                             --skip-convergence \
                             --cpus 8
    '''

}

process construct_boonstim_outputs{

    input:
        tuple val(sub), \
        path(msh), path(t1fs), path(m2m), path(fs), \
        path(l_pial), path(r_pial), \
        path(l_white), path(r_white), \
        path(l_thick), path(r_thick), \
        path(l_msm), path(r_msm), \
        path("${sub}.weightfunc.dscalar.nii"), path("${sub}.mask.dscalar.nii"), \
        path(C), path(R), path(bounds), path(qc_param), \
        path(femfunc)
//        path(loc), path(rot), path(hist)

    output:
        tuple val(sub), path("$sub"), emit: subject

    shell:
    '''
    #!/bin/bash

    # Make subdirectories
    mkdir !{sub}
    mkdir surfaces
    mkdir optimization

    # Move surface files
    mv \
        !{l_pial} !{r_pial} !{l_white} !{r_white} !{l_thick} !{r_thick} \
        !{l_msm} !{r_msm} !{sub}.weightfunc.dscalar.nii !{sub}.mask.dscalar.nii surfaces

    # Move optimization files
    #mv {loc} {rot} {hist} optimization
    mv !{C} !{R} !{bounds} !{femfunc} optimization

    # Move all files
    mv * !{sub} || true
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0],it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
            }

workflow {

    main:

        // Main preprocessing routine
        cifti_mesh_result = cifti_meshing(bids_channel)
        make_giftis_result = make_giftis(cifti_mesh_result.mesh_fs)
        registration_wf(cifti_mesh_result.mesh_fs)

        // User-defined weightfunction workflow
        weightfunc_input = cifti_mesh_result.fmriprep
                                            .join( cifti_mesh_result.cifti, by : 0 )
        weightfunc_input
        weightfunc_wf(weightfunc_input)

        // Calculate centroid on resampled data
        centroid_mask_input = weightfunc_wf.out.weightfunc
                                            .join(weightfunc_wf.out.weightfunc)
        centroid_mask(centroid_mask_input)
        resamplemask_wf(centroid_mask.out.masked, registration_wf.out.msm_sphere)
        centroid_wf(resamplemask_wf.out.resampled,
                    make_giftis_result.pial,
                    make_giftis_result.white,
                    make_giftis_result.midthickness,
                    cifti_mesh_result.t1fs_conform)

        // Parameterize the surface
        parameterization_wf(cifti_mesh_result.msh, centroid_wf.out.centroid)

        // Tetrahedral workflow
        dilate_mask_input = weightfunc_wf.out.weightfunc
                                            .join(cifti_mesh_result.cifti, by: 0)
                                            .map{ s,w,c ->  [
                                                                s,w,
                                                                "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                                                "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.midthickness.32k_fs_LR.surf.gii"
                                                            ]
                                                }
        dilate_mask(dilate_mask_input)
        weightfunc_mask_input = weightfunc_wf.out.weightfunc
                                            .join(dilate_mask.out.dilated, by: 0)
        weightfunc_mask(weightfunc_mask_input)
        resampleweightfunc_wf(weightfunc_mask.out.masked, registration_wf.out.msm_sphere)
        tet_project_wf(resampleweightfunc_wf.out.resampled,
                        make_giftis_result.pial,
                        make_giftis_result.white,
                        make_giftis_result.midthickness,
                        cifti_mesh_result.t1fs_conform,
                        cifti_mesh_result.msh)

        //// Gather inputs for optimization
        //optimize_inputs = cifti_mesh_result.msh
        //                            .join(tet_project_wf.out.fem_weights, by: 0)
        //                            .join(parameterization_wf.out.C, by: 0)
        //                            .join(parameterization_wf.out.R, by: 0)
        //                            .join(parameterization_wf.out.bounds, by: 0)
        //                            .map{ s,m,w,C,R,b -> [ s,m,w,C,R,b,"$params.coil" ] }
        //optimize_coil(optimize_inputs)


        //// Set up outputs
        registration_wf.out.msm_sphere.branch(lr_branch).set { msm }
        make_giftis_result.pial.branch(lr_branch).set { pial }
        make_giftis_result.white.branch(lr_branch).set { white }
        make_giftis_result.midthickness.branch(lr_branch).set { midthick }
        construct_output_input = cifti_mesh_result.msh
                                    .join(cifti_mesh_result.t1fs_conform, by: 0)
                                    .join(cifti_mesh_result.mesh_m2m, by: 0)
                                    .join(cifti_mesh_result.mesh_fs, by: 0)
                                    .join(pial.left, by:0)
                                    .join(pial.right, by:0)
                                    .join(white.left, by:0)
                                    .join(white.right, by:0)
                                    .join(midthick.left, by:0)
                                    .join(midthick.right, by:0)
                                    .join(msm.left, by:0)
                                    .join(msm.right, by:0)
                                    .join(resampleweightfunc_wf.out.resampled, by:0)
                                    .join(resamplemask_wf.out.resampled)
                                    .join(parameterization_wf.out.C)
                                    .join(parameterization_wf.out.R)
                                    .join(parameterization_wf.out.bounds)
                                    .join(parameterization_wf.out.qc_param)
                                    .join(tet_project_wf.out.fem_weights)
        //                            .join(weightfunc_wf.out.mask, by: 0)
        //                            .join(optimize_coil.out.position, by: 0)
        //                            .join(optimize_coil.out.orientation, by: 0)
        //                            .join(optimize_coil.out.history, by: 0)
        construct_boonstim_outputs(construct_output_input)

        publish:
            cifti_mesh_result.cifti   to: "$params.out/ciftify", mode: 'copy'
            cifti_mesh_result.fmriprep to: "$params.out/fmriprep", mode: 'copy'
            cifti_mesh_result.freesurfer to: "$params.out/freesurfer", mode: 'copy'
            cifti_mesh_result.fmriprep_html to: "$params.out/fmriprep", mode: 'copy'
            construct_boonstim_outputs.out.subject to: "$params.out/boonstim", mode: 'copy'

}

// TODO BUILD QC OUTPUTS
// TODO INSPECT AND MODIFY INPUT FILES IF NEEDED?
// TODO BUILD INPUT SUBSTITUTEABLE VERSION
// TODO LOGGING

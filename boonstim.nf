nextflow.preview.dsl=2

usage = file("${workflow.scriptFile.getParent()}/usage")
bindings = ["subjects": "$params.subjects",
            "cache_dir": "$params.cache_dir",
            "fmriprep": "$params.fmriprep",
            "ciftify": "$params.ciftify",
            "connectome": "$params.connectome",
            "bin": "$params.bin",
            "coil": "$params.coil",
            "license": "$params.license",
            "fmriprep_invocation": "$params.fmriprep_invocation",
            "fmriprep_anat_invocation": "$params.fmriprep_anat_invocation",
            "fmriprep_descriptor": "$params.fmriprep_descriptor",
            "ciftify_invocation": "$params.ciftify_invocation",
            "ciftify_descriptor": "$params.ciftify_descriptor",
            "weightworkflow" : "$params.weightworkflow"]
engine = new groovy.text.SimpleTemplateEngine()
toprint = engine.createTemplate(usage.text).make(bindings)
printhelp = params.help

req_param = ["--bids": "$params.bids",
             "--out": "$params.out"]
req_config_param = [
                    "fmriprep": "$params.fmriprep",
                    "ciftify": "$params.ciftify",
                    "connectome": "$params.connectome",
                    "bin": "$params.bin",
                    "coil": "$params.coil",
                    "license": "$params.license",
                    "fmriprep_invocation": "$params.fmriprep_invocation",
                    "fmriprep_anat_invocation": "$params.fmriprep_anat_invocation",
                    "fmriprep_descriptor": "$params.fmriprep_descriptor",
                    "ciftify_invocation": "$params.ciftify_invocation",
                    "ciftify_descriptor": "$params.ciftify_descriptor",
                    "weightworkflow" : "$params.weightworkflow"
                   ]
missing_req = req_param.grep{ (it.value == null || it.value == "") }
missing_req_config = req_config_param.grep{ (it.value == null || it.value == "") }

if (missing_req){
    log.error("Missing required command-line argument(s)!")
    missing_req.each{ log.error("Missing ${it.key}") }
    printhelp = true
}

if (missing_req_config){
    log.error("Config file missing required parameter(s)!")
    missing_req_config.each{ log.error("Please fill ${it.key} in config") }
    printhelp = true
}

if (printhelp){
    print(toprint)
    System.exit(0)
}

log.info("BIDS Directory: $params.bids")
log.info("Output Directory: $params.out")
if (params.subjects) {
    log.info ("Subject list file provided: $params.subjects")
}

log.info("Using Descriptor Files: $params.fmriprep_descriptor and $params.ciftify_descriptor")
log.info("Using Invocation Files: $params.fmriprep_invocation, $params.fmriprep_anat_invocation and $params.ciftify_invocation")
log.info("Using containers: $params.fmriprep and $params.ciftify")
log.info("Using user-defined ROI workflow: $params.weightworkflow")

///////////////////////////////////////////////////////////////////////////////

// IMPORT MODULES WORKFLOWS
include {cifti_meshing_wf as cifti_mesh_wf} from './modules/cifti_mesh_wf.nf' params(params)
include {make_giftis} from './modules/fs2gifti.nf' params(params)
include {registration_wf} from './modules/register_fs2cifti_wf.nf' params(params)
include {weightfunc_wf} from "${params.weightworkflow}" params(params)
include {resample2native_wf as resamplemask_wf} from './modules/resample2native.nf' params(params)
include {resample2native_wf as resampleweightfunc_wf} from './modules/resample2native.nf' params(params)
include {resample2native_wf as resampledistmap_wf} from './modules/resample2native.nf' params(params)
include {resample2native_wf as resamplesulc_wf} from './modules/resample2native.nf' params(params)
include {centroid_radial_wf as centroid_wf} from './modules/centroid_wf.nf' params(params)
include {tet_project_wf as tet_project_weightfunc_wf} from './modules/tetrahedral_wf.nf' params(params)
include {tet_project_wf as tet_project_roi_wf} from './modules/tetrahedral_wf.nf' params(params)
include {calculate_reference_field_wf} from './modules/reference_field_wf.nf' params(params)
include {fieldscaling_wf} from './modules/field_scaling.nf' params(params)
include {optimize_wf} from "./modules/optimization.nf" params(params)
include {qc_wf} from "./modules/qc.nf" params(params)

// IMPORT MODULES PROCESSES
include {apply_mask as centroid_mask} from './modules/utils.nf' params(params)
include {apply_mask as weightfunc_mask} from './modules/utils.nf' params(params)
include {cifti_dilate as dilate_mask} from './modules/utils.nf' params(params)

///////////////////////////////////////////////////////////////////////////////

//// Extract subjects to run
all_dirs = file(params.bids).list()
input_dirs = new File(params.bids).list()
output_dirs = new File(params.out).list()

input_channel = Channel.fromPath("$params.bids/sub-*", type: 'dir')
                    .map{i -> i.getBaseName()}

if (params.subjects){
    subjects_channel = Channel.fromPath(params.subjects)
                            .splitText(){it.strip()}
    input_channel = input_channel.join(subjects_channel)
}

if (!params.rewrite){
    out_channel = Channel.fromPath("$params.out/boonstim/sub-*", type: 'dir')
                    .map{o -> [o.getBaseName(), "o"]}
                    .ifEmpty(["", "o"])

    input_channel = input_channel.join(out_channel, remainder: true)
                        .filter{it.last() == null}
                        .map{i,n -> i}
}

process publish_base{

    publishDir path: "${params.out}/boonstim/${sub}", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(t1fs), path(centroid),\
    path(fem), path(dscalar)

    output:
    tuple path(t1fs), path(centroid),\
    path(fem), path(dscalar),\
    emit: base_out

    shell:
    '''
    #!/bin/bash

    echo "Moving files into boonstim/!{sub}..."
    '''

}

process publish_surfs{

    publishDir path: "${params.out}/boonstim/${sub}/T1w", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(pl), path(pr),\
    path(wl), path(wr),\
    path(ml), path(mr),\
    path(msml), path(msmr)

    output:
    tuple path(pl), path(pr),\
    path(wl), path(wr),\
    path(ml), path(mr),\
    path(msml), path(msmr), emit: surfs_out

    shell:
    '''
    #!/bin/bash
    echo "Transferring surfaces to boonstim/!{sub}/T1w..."
    '''

}

process publish_mri2mesh{

    publishDir path: "${params.out}/boonstim/${sub}", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(m2m), path(fs)

    output:
    tuple val(sub), path(m2m), path(fs)

    shell:
    '''
    #!/bin/bash
    echo "Transferring !{m2m} and !{fs} to boonstim/!{sub} folder..."
    '''
}

process publish_opt{

    publishDir path: "${params.out}/boonstim/${sub}/results", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(fields), path(coil), path(history),\
    path(brainsight), path(localite)

    output:
    tuple path(fields), path(coil), path(history),\
    path(brainsight), path(localite)

    shell:
    '''
    #!/bin/bash
    echo "Transferring optimization results to boonstim/!{sub}/results..."
    '''
}

process publish_scaleref{

    publishDir path: "${params.out}/boonstim/${sub}", \
               pattern: "*scalefactor.txt", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    val(name), path(factor)

    output:
    tuple val(sub), path("${sub}.${name}_scalefactor.txt")

    shell:
    '''
    echo "Moving stimulation scaling factor values into boonstim/!{sub}..."
    cp !{factor} "!{sub}.!{name}_scalefactor.txt"
    '''

}

process publish_cifti{

    stageInMode 'copy'
    publishDir path: "$params.out", \
               mode: 'move'

    input:
    tuple val(sub), \
    path("ciftify/*"), path("fmriprep/*"), path("freesurfer/*"), \
    path("ciftify/zz_templates")

    output:
    tuple path("ciftify/${sub}"), path("fmriprep/${sub}"),\
    path("freesurfer/${sub}")

    shell:
    '''
    echo "Copying fMRIPrep_Ciftify outputs"
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0], it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
                }

workflow {

    main:

        // Main preprocessing routine
        cifti_mesh_wf(input_channel)
        make_giftis(cifti_mesh_wf.out.mesh_fs)
        registration_wf(cifti_mesh_wf.out.mesh_fs)

        // User-defined weightfunction workflow
        weightfunc_input = cifti_mesh_wf.out.fmriprep
                                            .join( cifti_mesh_wf.out.cifti, by : 0 )
        weightfunc_input
        weightfunc_wf(weightfunc_input)

        // Calculate centroid on resampled data
        centroid_mask_input = weightfunc_wf.out.weightfunc
                                            .join(weightfunc_wf.out.mask)
        centroid_mask(centroid_mask_input)
        resamplemask_wf(centroid_mask.out.masked, registration_wf.out.msm_sphere)

        // Tetrahedral workflow

        // Resample the weightfunction
        dilate_mask_input = weightfunc_wf.out.mask
                                            .join(cifti_mesh_wf.out.cifti, by: 0)
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

        // Calculate a scalp seed
        centroid_wf(cifti_mesh_wf.out.msh,
                    resampleweightfunc_wf.out.resampled,
                    make_giftis.out.pial)

        tet_project_weightfunc_wf(resampleweightfunc_wf.out.resampled,
                        make_giftis.out.pial,
                        make_giftis.out.white,
                        make_giftis.out.midthickness,
                        cifti_mesh_wf.out.t1fs_conform,
                        cifti_mesh_wf.out.msh)

        // Gather inputs for optimization
        optimize_wf(
                    cifti_mesh_wf.out.msh,
                    tet_project_weightfunc_wf.out.fem_weights,
                    centroid_wf.out.centroid,
                    params.coil
                   )

        // Calculate scaling factor between coil and cortex across multiple references
        calculate_reference_field_wf(cifti_mesh_wf.out.cifti,
                                     Channel.from(params.ref_coords))
        resampledistmap_wf(
            calculate_reference_field_wf.out.rois,
            registration_wf.out.msm_sphere
        )

        fieldscaling_wf(
                        optimize_wf.out.fields,
                        make_giftis.out.pial,
                        resampledistmap_wf.out.resampled,
                        optimize_wf.out.matsimnibs
                      )

        // Gather BOONStim outputs for publishing
        registration_wf.out.msm_sphere.branch(lr_branch).set { msm }
        make_giftis.out.pial.branch(lr_branch).set { pial }
        make_giftis.out.white.branch(lr_branch).set { white }
        make_giftis.out.midthickness.branch(lr_branch).set { midthick }

        // Run surface QC workflow
        i_resample_sulc = cifti_mesh_wf.out.cifti
                            .map{s,c -> [s,
                                        "${c}/MNINonLinear/fsaverage_LR32k/${s}.sulc.32k_fs_LR.dscalar.nii"
                                        ]
                                }

        resamplesulc_wf(
            i_resample_sulc,
            registration_wf.out.msm_sphere
        )


        // Generate QC images
        //qc_wf(
        //    resampleweightfunc_wf.out.resampled,
        //    pial.left, pial.right,
        //    resamplesulc_wf.out.resampled)

        /* Step 1: Publish base outputs */
        i_publish_base = cifti_mesh_wf.out.t1fs_conform
                                .join(centroid_wf.out.centroid)
                                .join(tet_project_weightfunc_wf.out.fem_weights)
                                .join(resampleweightfunc_wf.out.resampled)
        publish_base(i_publish_base)

        /* Step 2: Publish native space surfaces used to map out
        weight function
        */
        i_publish_surfs = pial.left.join(pial.right)
                            .join(white.left).join(white.right)
                            .join(midthick.left).join(midthick.right)
                            .join(msm.left).join(msm.right)
        publish_surfs(i_publish_surfs)

        /* Step 3: Publish meshing results from mri2mesh */
        i_publish_mri2mesh = cifti_mesh_wf.out.mesh_m2m
                                .join(cifti_mesh_wf.out.mesh_fs)
        publish_mri2mesh(i_publish_mri2mesh)

        /* Step 4: Publish optimization results */
        i_publish_opt = optimize_wf.out.fields
                            .join(optimize_wf.out.coil)
                            .join(optimize_wf.out.history)
                            .join(optimize_wf.out.brainsight)
                            .join(optimize_wf.out.localite)
        publish_opt(i_publish_opt)

        /* Step 5: Publish the reference scaling values */
        publish_scaleref(fieldscaling_wf.out.scaling_factor)

        // Publish Ciftify outputs
        publish_cifti_input = cifti_mesh_wf.out.cifti
                                    .join(cifti_mesh_wf.out.fmriprep)
                                    .join(cifti_mesh_wf.out.freesurfer)
                                    .combine(["$params.zz"])
        publish_cifti(publish_cifti_input)
}

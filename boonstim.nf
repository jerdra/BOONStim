nextflow.preview.dsl=2

include { getArgumentParser } from "./lib/args"

parser = getArgumentParser(
    title: "BOONStim TMS Optimization Pipeline",
    description: "An end-to-end pipeline for integrating fMRI data into TMS target\
 derivation and optimization",
    scriptName: "${workflow.scriptName}".toString(),
    note: "Configuration arguments should be defined in a .nf.config file (use -c arg), or\
 a .json file that can be used (use -params-file arg)"
)

parser.addRequired("--bids",
    "Path to BIDS directory",
    params.bids.toString(),
    "BIDS_DIRECTORY")

parser.addRequired("--out",
    "Path to output directory",
    params.out.toString(),
    "OUTPUT_DIR")

parser.addConfigOpt("--connectome",
    "Path to Connectome Workbench Image",
    params.connectome.toString(),
    "CONNECTOME_IMG")

parser.addConfigOpt("--bin",
    "Path to BOONStim bin directory",
    params.bin.toString(),
    "BIN_DIR")

parser.addConfigOpt("--coil",
    "Path to COIL file (.ccd or .nii.gz)",
    params.coil.toString(),
    "COIL_FILE")

parser.addConfigOpt("--license",
    "Path to Freesurfer license file",
    params.license.toString(),
    "FS_LICENSE")

parser.addConfigOpt("--fmriprep",
    "Path to fMRIPrep Image",
    params.fmriprep.toString(),
    "FMRIPREP_IMG")

parser.addConfigOpt("--fmriprep_descriptor",
    "Path to fMRIPrep Descriptor File",
    params.fmriprep_descriptor.toString(),
    "FMRIPREP_DESCRIPTOR")

parser.addConfigOpt("--fmriprep_invocation",
    "Path to fMRIPrep invocation file",
    params.fmriprep_invocation.toString(),
    "FMRIPREP_INVOCATION")

parser.addConfigOpt("--anat_invocation",
    "Path to fMRIPrep anatomical only invocation file",
    params.anat_invocation.toString(),
    "FMRIPREP_ANAT_INVOCATION")

parser.addConfigOpt("--ciftify",
    "Path to Ciftify Image",
    params.ciftify.toString(),
    "CIFTIFY_IMG")

parser.addConfigOpt("--ciftify_invocation",
    "Path to Ciftify invocation file",
    params.ciftify_invocation.toString(),
    "CIFTIFY_INVOCATION")

parser.addConfigOpt("--ciftify_descriptor",
    "Path to Ciftify descriptor file",
    params.ciftify_descriptor.toString(),
    "CIFTIFY_DESCRIPTOR")

parser.addConfigOpt("--simnibs",
    "Path to SimNIBS container image",
    params.simnibs.toString(),
    "SIMNIBS_IMG")

parser.addConfigOpt("--fieldopt",
    "Path to FieldOpt container image",
    params.fieldopt.toString(),
    "FIELDOPT_IMG")

parser.addConfigOpt("--weightworkflow",
    "Path to weight workflow module entrypoint",
    params.weightworkflow.toString(),
    "WEIGHTWORKFLOW_ENTRYPOINT")

parser.addOptional("--subjects",
    "Path to subject text file containing 1 BIDS subject/line",
    "SUBJECT_FILE")

parser.addOptional("--coil-distances",
    "Path to text file containing [subject, coil distance] tuples"
    "COIL_DISTANCES")

parser.addOptional("--num_cpus",
    "Maximum number of threads to use when submitting jobs [Default: $params.num_cpus]",
    "NUM_CPUS")

parser.addOptional("--cache_dir",
    "Create a cache directory to store intermediate results to speed up reruns",
    "CACHE_DIR")

missingArgs = parser.isMissingRequired()
missingConfig = parser.isMissingConfig()

if (params.help) {
    print(parser.makeDoc())
    System.exit(0)
}

if (missingArgs || missingConfig) {
    log.error("Missing required parameters")
    missingArgs.each{ log.error("Missing ${it}") }
    missingConfig.each{ log.error("Missing ${it}") }
    print(parser.makeDoc())
    System.exit(1)
}

include {weightfunc_wf} from "${params.weightworkflow}" params(params)
include {cifti_meshing_wf as cifti_mesh_wf} from './modules/cifti_mesh_wf.nf' params(params)
include {make_giftis} from './modules/fs2gifti.nf' params(params)
include {registration_wf} from './modules/register_fs2cifti_wf.nf' params(params)
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
include {apply_mask as centroid_mask} from './modules/utils.nf' params(params)
include {apply_mask as weightfunc_mask} from './modules/utils.nf' params(params)
include {cifti_dilate as dilate_mask} from './modules/utils.nf' params(params)

log.info("BIDS Directory: $params.bids")
log.info("Output Directory: $params.out")
if (params.subjects) {
    log.info ("Subject list file provided: $params.subjects")
}

log.info("Using Descriptor Files: $params.fmriprep_descriptor and $params.ciftify_descriptor")
log.info("Using Invocation Files: $params.fmriprep_invocation, $params.fmriprep_anat_invocation and $params.ciftify_invocation")
log.info("Using containers: $params.fmriprep and $params.ciftify")
log.info("Using user-defined ROI workflow: $params.weightworkflow")

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

    publishDir path: "${params.out}/boonstim/${sub}/results", \
               pattern: "${sub}.${name}*", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    val(name), path(factor), path(html), path(geo)

    output:
    tuple val(sub), path("${sub}.${name}*")

    shell:
    '''
    echo "Moving stimulation scaling factor values into boonstim/!{sub}..."
    cp !{factor} "!{sub}.!{name}_scalefactor.txt"
    cp !{html} "!{sub}.!{name}_qc.html"
    cp !{geo} "!{sub}.!{name}_qc.geo"
    '''

}

process publish_cifti{

    stageInMode 'copy'

    publishDir path: "$params.out",\
               mode: 'move'

    input:
    tuple val(sub),\
    path("ciftify/${sub}"), path("ciftify/qc_fmri/*"),\
    path("ciftify/qc_recon_all/${sub}"),\
    path("fmriprep/*"), path(html),\
    path("freesurfer/*"),\
    path("ciftify/zz_templates")

    output:
    tuple path("ciftify/${sub}"),
    path("ciftify/qc_fmri/${sub}*", includeInputs: true),\
    path("ciftify/qc_recon_all/${sub}"),\
    path("fmriprep/${sub}"), path("fmriprep/${sub}.html"),\
    path("freesurfer/${sub}"), path("ciftify/zz_templates")

    shell:
    '''
    echo "Copying fMRIPrep_Ciftify outputs"
    mv !{html} fmriprep/

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
        i_publish_scaleref = fieldscaling_wf.out.scaling_factor
                                            .join(fieldscaling_wf.out.qc_html, by: [0,1])
                                            .join(fieldscaling_wf.out.qc_geo, by: [0,1])
        publish_scaleref(i_publish_scaleref)

        // Publish Ciftify outputs
        publish_cifti_input = cifti_mesh_wf.out.cifti
                                    .join(cifti_mesh_wf.out.cifti_qc_fmri)
                                    .join(cifti_mesh_wf.out.cifti_qc_recon)
                                    .join(cifti_mesh_wf.out.fmriprep)
                                    .join(cifti_mesh_wf.out.fmriprep_html)
                                    .join(cifti_mesh_wf.out.freesurfer)
                                    .combine(["$params.zz"])
        publish_cifti(publish_cifti_input)
}

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

include {scalar_optimization} from "./pipeline/scalar_optimization.nf" params(params)
include {coordinate_optimization} from "./pipeline/coordinate_optimization.nf" params(params)

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

workflow scalar_workflow {
    scalar_optimization(input_channel)
}


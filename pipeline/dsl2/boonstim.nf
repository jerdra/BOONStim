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

log.info("Using Descriptor Files: $params.anat_descriptor and $params.ciftify_descriptor")
log.info("Using Invocation Files: $params.anat_invocation and $params.ciftify_invocation")


///////////////////////////////////////////////////////////////////////////////

// IMPORT WORKFLOWS
include cifti_meshing from './modules/cifti_mesh_wf.nf' params(params)

///////////////////////////////////////////////////////////////////////////////

// Main Workflow

all_dirs = file(params.bids)
bids_channel = Channel
                    .from(all_dirs.list())
                    .filter { it.contains('sub') }
                    .take(1)

// Process definitions

workflow {

    cifti_meshing(bids_channel)

}

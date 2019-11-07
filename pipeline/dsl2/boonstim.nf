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
include make_giftis from './modules/fs2gifti.nf' params(params)
include registration_wf from './modules/register_fs2cifti_wf.nf' params(params)
include weightfunc_wf from "${params.weightworkflow}" params(params)

///////////////////////////////////////////////////////////////////////////////

// Main Workflow

all_dirs = file(params.bids)
bids_channel = Channel
                    .from(all_dirs.list())
                    .filter { it.contains('sub-CMH') }
                    .take(1)

// Process definitions

workflow {

    main:
        
        // Main preprocessing routine
        cifti_mesh_result = cifti_meshing(bids_channel)

        // Space registration
        make_giftis_result = make_giftis(cifti_mesh_result.mesh_fs)
        registration_wf(cifti_mesh_result.mesh_fs)

        // Calculation of custom weightfunction and centroid
        weightfunc_input = cifti_mesh_result.fmriprep
                                            .join( cifti_mesh_result.cifti, by : 0 )
        weightfunc_wf(weightfunc_input)

        // Resample the data into native space


        //Check outputs for debugging
        weightfunc_wf.out.weightfunc | view
        weightfunc_wf.out.mask | view
        make_giftis_result.pial | view
        make_giftis_result.white | view
        make_giftis_result.midthickness | view

        


}

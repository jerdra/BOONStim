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
include resample2native_wf as resamplemask_wf from './modules/resample2native.nf' params(params)
include resample2native_wf as resampleweightfunc_wf from './modules/resample2native.nf' params(params)
include centroid_wf from './modules/centroid_wf.nf' params(params)
include parameterization_wf from './modules/surfparams_wf.nf' params(params)
include tet_project_wf from './modules/tetrahedral_wf.nf' params(params)

///////////////////////////////////////////////////////////////////////////////

// Main Workflow

all_dirs = file(params.bids)
bids_channel = Channel
                    .from(all_dirs.list())
                    .filter { it.contains('sub-CMH') }
                    .take(1)

// Process definitions
process optimize_coil {

    label 'numpy'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh), path(weights), path(C), path(R), path(bounds), path(coil)

    output:
    tuple val(sub), path('coil_position'), path('coil_orientation')

    shell:
    ''' 
    /scripts/optimize_fem.py !{msh} !{weights} !{C} !{bounds} !{R} !{coil} \
                             coil_position coil_orientation
    '''
    
}

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

        // Resample the mask
        resamplemask_wf(weightfunc_wf.out.mask, registration_wf.out.msm_sphere)

        // Calculate centroid on resampled data
        centroid_wf(resamplemask_wf.out.resampled,make_giftis_result.pial, make_giftis_result.white, make_giftis_result.midthickness, cifti_mesh_result.t1fs_conform)

        // Parameterize the surface
        parameterization_wf(cifti_mesh_result.msh, centroid_wf.out.centroid) 
        // Resample the weightfunc
        resampleweightfunc_wf(weightfunc_wf.out.weightfunc, registration_wf.out.msm_sphere)
        
        // Project weightfunc into tetrahedral mesh space
        //tet_project_wf(resampleweightfunc_wf.out.resampled, make_giftis_result.pial, make_giftis_result.white, make_giftis_result.midthickness, cifti_mesh_result.t1fs_conform, cifti_mesh_result.msh)

        //// Gather inputs for optimization
        //optimize_inputs = cifti_mesh_result.msh
        //                            .join(tet_project_wf.out.fem_weights, by: 0)
        //                            .join(parameterization_wf.out.C, by: 0)
        //                            .join(parameterization_wf.out.R, by: 0)
        //                            .join(parameterization_wf.out.bounds, by: 0)
        //                            .map{ s,m,w,C,R,b -> [ s,m,w,C,R,b,"$params.coil" ] }
        //optimize_coil(optimize_inputs)

        // Formulate outputs for ciftify

        publish:
            
            // Ciftify pipeline outputs
            cifti_mesh_result.cifti   to: "$params.out", mode: 'copy', pattern: {"./ciftify/*"}
            cifti_mesh_result.fmriprep to: "$params.out", mode: 'copy', pattern: {"./fmriprep/*"}
            cifti_mesh_result.freesurfer to: "$params.out", mode: 'copy', pattern: {"./freesurfer/*"}
            cifti_mesh_result.freesurfer to: "$params.out", mode: 'copy', pattern: {"./freesurfer/${it[0]}"}
            // Meshing outptus
            

}

// TODO LIST:
// TODO BUILD QC OUTPUTS
// TODO INSPECT AND MODIFY INPUT FILES IF NEEDED?
// TODO PUBLISH OUTPUTS

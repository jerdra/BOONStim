nextflow.preview.dsl=2

include { weightfunc_wf } from "${params.weightworkflow}" params(params)
include {cifti_meshing_wf as cifti_mesh_wf} from '../modules/cifti_mesh_wf.nf' params(params)
include { rigidRegistration } from "../modules/utils.nf"

// TODO: Include ADM optimization module


// Default params 
params.radius = 20 // in mm

process coordinate_transform {

    /*
    * Order of transformations to points:
    * 1. Apply T1w --> MNI nonlinear transform no Affine
    * 2. Apply T1w --> MNI Affine transform
    * 3. Apply SimT1 --> CiftiT1 Rigid body transform
    * Alternatively can run surface-based resampling between
    * CiftiMNI and CiftiT1 space
    * Then just run RigidBody
    */

    // TODO: Test series of coordinate transforms on a sample ROI
    // TODO: Have Surface resampling as a back-up option

    input:
    tuple val(sub), path(coords), path(transforms)

    output:
    tuple val(sub), path("${sub}_target_coord.txt"), emit: transformed

    shell:
    """
    
    """
}

workflow coordinate_optimization {

    take:
    subject_channel

    main:
    cifti_meshing_wf(subject_channel)

    weightfunc_input = cifti_mesh_wf.out.fmriprep
                                        .join( cifti_mesh_wf.out.cifti, by : 0 )
    weightfunc_input
    weightfunc_wf(weightfunc_input)

    i_rigidRegistration = cifti_meshing_wf.out.cifti
        .join(cifti_meshing_wf.out.mesh_m2m, by: 0)
        .map{ s, c, m2m -> [
            s,
            "${c}/T1w/T1w.nii.gz",
            "${m2m}/T1fs_resamp.nii.gz"
            ]
        }

    rigidRegistration(i_rigidRegistration)

    i_coordinate_transform = cifti_meshing_wf.out.cifti
        .map { s, c -> [
            s,
            "${c}/MNINonLinear/xfms/T1w2StandardLinear.mat",
            "${c}/MNINonLinear/xfms/T1w2Standard_warp_noaffine.nii.gz"
        ]}
        .join(rigidRegistration.out.transform)
    coordinate_transform(i_coordinate_transform)

    // TODO: Feed coordinate into ADM pipeline
    i_adm = cifti_meshing_wf.out.msh
        .join(coordinate_transform.out.transformed)
        .spread([params.radius])
    // adm_wf(i_adm)
}

nextflow.preview.dsl=2
params.radius = 20 // in mm

include { weightfunc_wf } from "${params.weightworkflow}" params(params)
include { cifti_meshing_wf as cifti_mesh_wf } from '../modules/cifti_mesh_wf.nf' params(params)
include { rigidRegistration; coordinate_transform; map_coordinate } from "../modules/transformation.nf" params(params)

workflow coordinate_optimization {

    /*
    Perform single-coordinate optimization with a set radius

    Arguments:
        subject_channel (channel): (subject_id: str, subject_params: Map) Subject ID and additional parameters to be passed to optimization. 
            Supported parameters: ['hair_thickness']

    Parameters:
        radius (float): Radius of optimization target around derived coordinate
    */

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
    map_coordinate(weightfunc_wf.out.coordinate, rigidRegistration.out.transform)

    i_adm = cifti_meshing_wf.out.msh
        .join(cifti_meshing_wf.out.mesh_fs)
        .join(map_coordinate.out.ras_coords)
        .spread([params.radius])
    adm_wf(i_adm)
}

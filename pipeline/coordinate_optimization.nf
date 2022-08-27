nextflow.preview.dsl=2
params.radius = 20 // in mm

include { weightfunc_wf } from "${params.weightworkflow}" params(params)
include { cifti_meshing_wf as cifti_mesh_wf } from '../modules/cifti_mesh_wf.nf' params(params)
include { rigidRegistration; coordinate_transform; map_coordinate } from "../modules/transformation.nf" params(params)
include { dosage_adjustment_wf } from "../modules/dosage_wf.nf" params(params)
include { adm_wf } from "../modules/adm_opt.nf" params(params)
include { neuronav_wf } from "../modules/neuronav.nf" params(params)

def try_as_numeric(map){

    /*
    Attempt to transform values of a map
    to numeric values
    */

    map.collectEntries { key, value ->
        try {
            def new_value = value as float
        } catch(NumberFormatException e){
            def new_value = value
        }
        [key, value]
    }
}

workflow coordinate_optimization {

    /*
    Perform single-coordinate optimization with a set radius

    Arguments:
        subject_spec (channel): (subject_id: str, subject_params: Map) Subject ID and additional parameters to be passed to optimization. 
            Supported parameters: ['hair_thickness']

    Parameters:
        radius (float): Radius of optimization target around derived coordinate
    */

    take:
    subject_spec

    main:
    formatted_subject_spec = subject_spec.map{s, p -> [s, try_as_numeric(p)]}
    subject_channel = subject_spec.map { s, p -> s }
    cifti_mesh_wf(subject_channel)

    weightfunc_input = cifti_mesh_wf.out.fmriprep
                                         .join( cifti_mesh_wf.out.cifti, by : 0 )
    weightfunc_input
    weightfunc_wf(weightfunc_input)

    i_rigidRegistration = cifti_mesh_wf.out.cifti
        .join(cifti_mesh_wf.out.mesh_m2m, by: 0)
        .map{ s, c, m2m -> [
            s,
            "${c}/T1w/T1w.nii.gz",
            "${m2m}/T1fs_resamp.nii.gz"
            ]
        }

    rigidRegistration(i_rigidRegistration)
    map_coordinate(weightfunc_wf.out.coordinate, rigidRegistration.out.transform)
    adm_wf(
        formatted_subject_spec,
        cifti_mesh_wf.out.msh,
        cifti_mesh_wf.out.mesh_fs,
        map_coordinate.out.transformed,
        Channel.from(params.radius),
        Channel.fromPath(params.coil)
    )

    dosage_adjustment_wf(
        adm_wf.out.sim_msh,
        cifti_mesh_wf.out.mesh_m2m,
        params.reference_magnitude
    )

    neuronav_wf(adm_wf.out.matsimnibs)
}

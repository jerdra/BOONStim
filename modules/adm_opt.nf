nextflow.preview.dsl=2

import groovy.json.JsonOutput


// DEFAULT VARIABLE
params.optimize_magnitude = true

include { target_direction_wf } from '../modules/target_direction.nf' params(params)
include { numpy2txt } from "./utils.nf" params(params)


process prepare_parameters {
    /*
    Build DB storing parameters to be used for
    optimization

    Arguments:
        subject (str): Subject ID
        opt_spec (Map<str, Union[str,int,float,bool]>): Optimization specs
        target_spec (Map<str, float>): Target specification. The following
            keys are supported:
            ['pos_x', 'pos_y', 'pos_z', 'dir_x', 'dir_y', 'dir_z']

    Output:
        json (val): (json: Path) JSON containing parameters to use for optimization
    */

    input:
    tuple val(subject), val(settings)

    output:
    path("${subject}.json"), emit: json

    exec:


    // Add history empty map
    def final_map = settings + [history: [:]]
    def json_str = JsonOutput.toJson(final_map)
    def json_beauty = JsonOutput.prettyPrint(json_str)
    File file = new File(subject)
    file.write(json_beauty)
}

process adm_optimize {

    label 'fieldopt'

    /*
    Direction optimization variant of ADM

    Arguments:
        sub (str): Subject ID
        json (path): JSON containing inputs to ADM
        msh (Path): Path to .msh file
        radius (float): Optimization target radius

    Parameters:
        optimize_magnitude (bool): Whether to optimize for magnitude or not

    Outputs:
        sim_msh (channel): (sub, sim_msh: Path) Optimized E-field simulation .msh
        geo (channel): (sub, geo: Path) Optimized E-Field coil placement
        matsimnibs (channel): (sub, msn: Path) Optimal coil orientation matrix
    */

    input:
    tuple val(sub), path(json), path(msh), val(radius)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: sim_msh
    tuple val(sub), path("${sub}_optimized_fields.geo"), emit: geo
    tuple val(sub), path("${sub}_optimized_orientation.npy"), emit: matsimnibs

    script:
    
    def direction = (params.optimize_magnitude) ? "magnitude" : "direction"

    """
    /scripts/adm_optimize.py \
        ${msh} \
        ${json} \
        --radius ${radius} \
        --sim_result ${sub}_optimized_fields.msh \
        --sim_geo ${sub}_optimized_fields.geo \
        ${sub}_optimized_orientation.npy \
        ${direction}
    """
    
}



workflow adm_wf {

    /*
    Perform Auxiliary Dipole Method-based optimization on a single
    target coordinate

    Arguments:
      subject_spec (channel): (subject, params: Map) Subject parameters for
        running optimization
      msh (channel): (subject, mesh_file: Path)
      fs (channel): (subject, fs_dir: Path)
      coord (channel): (subject, coordinate: Path)
      radius (value): float

    Parameters:
      optimize_magnitude (bool): Whether to optimize the magnitude of the e-field
       or to optimize for the direction normal to the brain node closest to `coordinate`
       [default=true]

    Outputs:
        sim_msh (channel): (subject, sim_file: Path)
        geo (channel): (subject, geo_file: Path)
        matsimnibs (channel): (subject, msn_file: Path)

    */

    take:
    subject_spec
    msh
    fs
    coord
    radius

    main:

    // Prepare ADM parameters
    target_direction_wf(coord, fs)

    // Read coordinates into map
    coord_map = numpy2txt(coord)
        .map { sub, coords -> [
            sub,
            coords.text.strip("\n").split(",")
        ]}
        .map { sub, c -> [
            sub,
            [pos_x: c[0], pos_y: c[1], pos_z: c[2]]
        ]}

    direction_map = numpy2txt(target_direction_wf.out.direction)
        .map { sub, direction -> [
            sub,
            direction.text.strip("\n").split(",")
        ]}
        .map { sub, d -> [
            sub,
            [dir_x: d[0], dir_y: d[1], dir_z: d[2]]
        ]}

    // Merge target specs into single map
    target_spec = coord_map.join(direction_map)
        .map { s, c, d -> [s, c + d] }

    // Prepare ADM parameters into JSON file
    prepare_parameters(subject_spec.join(target_spec))
    adm_optimize( prepare_parameters.out.json.join(msh).spread([radius]) )

    emit:
    sim_msh = adm_optimization.out.sim_msh
    geo = adm_optimization.out.geo
    matsimnibs = adm_optimization.out.matsimnibs
}

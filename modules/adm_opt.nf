nextflow.preview.dsl=2

import groovy.json.JsonOutput


// DEFAULT VARIABLE
params.optimize_magnitude = true

include { target_direction_wf } from '../modules/target_direction.nf' params(params)
include { numpy2txt as coord2txt; numpy2txt as direction2txt } from "./utils.nf" params(params)


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
    tuple val(subject), val(settings), val(target_spec)

    output:
    tuple val(subject), path("${subject}.json"), emit: json

    exec:
    def final_map = settings + [history: []] + target_spec
    def json_str = JsonOutput.toJson(final_map)
    def json_beauty = JsonOutput.prettyPrint(json_str)
    File file = new File("${task.workDir}/${subject}.json")
    println task.workDir
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
    tuple val(sub), path(json), path(msh), val(radius), path(coil)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: sim_msh
    tuple val(sub), path("${sub}_optimized_fields.geo"), emit: geo
    tuple val(sub), path("${sub}_optimized_orientation.npy"), emit: matsimnibs

    script:

    def direction = (params.optimize_magnitude.toBoolean()) ? "magnitude" : "direction"

    """
    python /scripts/adm_optimize.py \
        ${msh} \
        ${json} \
        ${coil} \
        --radius ${radius} \
        --sim_result ${sub}_optimized_fields.msh \
        --sim_geo ${sub}_optimized_fields.geo \
        ${sub}_optimized_orientation.npy \
        ${direction}
    """

}

process make_adm_qc {
    /*
    Generate a quality control image from ADM outputs

    Arguments:
        sub (str): Subject ID
        sim_msh (Path): Simulation mesh containing fields + target
        msn (Path): Optimal coil position
        json (Path): Optimization settings containing coordinates and
            optionally, direction vectors

    Parameters:
        optimize_magnitude (bool): If optimize magnitude is enabled
            will omit showing direction vector as it would not
            have been used for optimization

    Outputs:
        qc_html (Path): Interactive QC HTML file
        qc_img (Path): Static QC image
    */
    label 'fieldopt'

    input:
    tuple val(sub), path(sim_msh), path(msn), path(json)

    output:
    tuple val(sub), path("${sub}_desc-qc.png"), emit: img
    tuple val(sub), path("${sub}_desc-qc.html"), emit: html

    script:
    def direction = (params.optimize_magnitude.toBoolean()) ? "magnitude" : "direction"

    """
    python /scripts/adm_qc.py \
        ${sim_msh} \
        ${msn} \
        ${json} \
        --export-img ${sub}_desc-qc.png \
        --export-html ${sub}_desc-qc.html \
        ${direction}
    """
}

process publish_adm {
    
    /*
    Publish outputs from ADM optimization and quality control
    into the output directory

    Arguments:
        sub (str): Subject ID
        sim_msh (Path): Simulation mesh
        sim_geo (Path): Coil dipole .geo
        matsimnibs (Path): Optimal coil placement matrix
        qc (Path): QC image for optimization
        spec (Path): Optimization settings used JSON file

    Parameters:
        out (Path): Base output directory
    */

    publishDir  path: "${params.out}/boonstim/${sub}/adm", \
                mode: 'copy', \
                overwrite: true

    input:
    tuple val(sub), path(sim_msh), path(sim_geo),\
    path(matsimnibs), path(qc_img), path(qc_html), path(spec),\
    path(dir_qc_img), path(dir_qc_html)

    output:
    tuple val(sub), path(sim_msh), path(sim_geo),\
    path(matsimnibs), path(qc_img), path(qc_html), path(spec),\
    path(dir_qc_img), path(dir_qc_html)

    shell:
    """
    echo 'Moving data into ${params.out}/boonstim/${sub}/adm!'
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
      coil (value): Coil definition file (.ccd) to use

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
    coil

    main:

    // Prepare ADM parameters
    target_direction_wf(coord, fs)

    // Read coordinates into map
    coord_map = coord2txt(coord)
        .map { sub, coords -> [
            sub,
            coords.text.split("\n")
        ]}
        .map { sub, c -> [
            sub,
            [pos_x: c[0], pos_y: c[1], pos_z: c[2]]
        ]}

   direction_map = direction2txt(target_direction_wf.out.direction)
       .map { sub, direction -> [
           sub,
           direction.text.split("\n")
       ]}
       .map { sub, d -> [
           sub,
           [dir_x: d[0], dir_y: d[1], dir_z: d[2]]
       ]}

   // Merge target specs into single map
   target_spec = coord_map.join(direction_map)
       .map { s, c, d -> [s, c + d] }
       .map { s, spec -> [
            s,
            spec.collectEntries{ key, value -> [key, value as float]}
        ]}

   // Prepare ADM parameters into JSON file
   prepare_parameters(subject_spec.join(target_spec))

   adm_optimize(
       prepare_parameters.out.json
           .join(msh)
           .combine(radius)
           .combine(coil)
   )

   // Generate a quality control image of optimization
   make_adm_qc(
    adm_optimize.out.sim_msh
        .join(adm_optimize.out.matsimnibs)
        .join(prepare_parameters.out.json)
   )

   publish_adm(
    adm_optimize.out.sim_msh
        .join(adm_optimize.out.geo)
        .join(adm_optimize.out.matsimnibs)
        .join(make_adm_qc.out.img)
        .join(make_adm_qc.out.html)
        .join(prepare_parameters.out.json)
        .join(target_direction_wf.out.qc_img)
        .join(target_direction_wf.out.qc_html)
   )


   emit:
   sim_msh = adm_optimize.out.sim_msh
   geo = adm_optimize.out.geo
   matsimnibs = adm_optimize.out.matsimnibs
   parameters = prepare_parameters.out.json
   qc_img = make_adm_qc.out.img
   qc_html = make_adm_qc.out.html
   direction_qc_img = target_direction_wf.out.qc_img
   direction_qc_html = target_direction_wf.out.qc_html
}

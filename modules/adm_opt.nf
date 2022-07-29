nextflow.preview.dsl=2

// DEFAULT VARIABLE
params.optimize_magnitude = true

include { target_direction_wf } from '../modules/target_direction.nf' params(params)


process adm_optimization_radial {

    label 'fieldopt'
    label 'bin'

    input:
    tuple val(sub), path(coord), path(direction), path(msh), val(radius)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: sim_msh
    tuple val(sub), path("${sub}_optimized_fields.geo"), emit: geo
    tuple val(sub), path("${sub}_optimized_orientation.npy"), emit: matsimnibs

    shell:
    """
    /scripts/adm_optimize.py \
        ${msh} \
        ${coord} \
        --direction ${direction} \
        --radius ${radius} \
        --sim_result ${sub}_optimized_fields.msh \
        --sim_geo ${sub}_optimized_fields.geo \
        ${sub}_optimized_orientation.npy
    """
}

process adm_optimization_mag {
    label 'fieldopt'
    label 'bin'

    input:
    tuple val(sub), path(coord), path(msh), val(radius)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: sim_msh
    tuple val(sub), path("${sub}_optimized_fields.geo"), emit: geo
    tuple val(sub), path("${sub}_optimized_orientation.npy"), emit: matsimnibs

    shell:
    """
    /scripts/adm_optimize.py \
        ${msh} \
        ${coord} \
        --radius ${radius} \
        --sim_result ${sub}_optimized_fields.msh \
        --sim_geo ${sub}_optimized_fields.geo \
        ${sub}_optimized_orientation.npy
    """
}


workflow adm_wf {

    /*
    Perform Auxiliary Dipole Method-based optimization on a single
    target coordinate
    
    Arguments:
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
    msh
    fs
    coord
    radius

    main:
    if (params.optimize_magnitude) {

        adm_optimization_mag(
            msh.join(coord) .spread([radius])
        )

    } else {

        target_direction_wf(coord, fs)
        adm_optimization_radial(
            msh.join(coord)
                .join(target_direction_wf.out.target_direction)
                .spread([radius])
        )

    }


    emit:
    sim_msh = adm_optimization.out.sim_msh
    geo = adm_optimization.out.geo
    matsimnibs = adm_optimization.out.matsimnibs
}

nextflow.preview.dsl=2

params.optimize_magnitude = true

process get_coordinate_normal {
    label 'fieldopt'
    label 'bin'

    input:
    tuple val(sub), path(coord), path(msh), path(fs)

    output:
    tuple val(sub), path("${sub}_target_normal.npy"), emit: target_direction

    shell:
    """
    /scripts/get_target_direction.py \
        ${fs} \
        ${msh} \
        ${coord} \
        ${sub}_target_normal.npy
    """
}

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

            get_coordinate_normal(coord.join(msh).join(fs)) 
            adm_optimization_radial(
                msh.join(coord)
                    .join(get_coordinate_normal.out.target_direction)
                    .spread([radius])
            )

        }


    emit:
        sim_msh = adm_optimization.out.sim_msh
        geo = adm_optimization.out.geo
        matsimnibs = adm_optimization.out.matsimnibs
}

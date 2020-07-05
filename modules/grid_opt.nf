nextflow.preview.dsl = 2

// Set up grid optimization schema
process grid_optimization{

    input:
    tuple val(sub), path(msh),\
    path(weights), path(centroid), path(coil)

    output:
    tuple val(sub), path("${sub}_orientation.txt"), emit: orientation
    tuple val(sub), path("${sub}_optimal_sim.msh"), emit: opt_msh
    tuple val(sub), path("${sub}_optimal_coilpos.geo"), emit: opt_coil
    tuple val(sub), path("${sub}_history.txt"), emit: history

    shell:
    '''
    #!/bin/bash

    #
    '''

}

workflow optimize{

    take:
       msh,
       weights,
       centroid,
       coil

    main:
        i_grid_optimization = msh.join(weights)
                                .join(centroid)
                                .join(coil)
        grid_optimization(i_grid_optimization)

    emit:
        orientation = optimize_coil.out.orientation
        opt_msh = optimize_coil.out.opt_msh
        opt_coil = optimize_coil.out.opt_coil
        history = optimize_coil.out.history
}

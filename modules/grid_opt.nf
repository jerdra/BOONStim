nextflow.preview.dsl = 2

// Set up grid optimization schema
process grid_optimization{

    stageInMode 'copy'
    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh),\
    path(weights), path(centroid), path(coil)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: fields
    tuple val(sub), path("${sub}_optimized_coil.geo"), emit: coil
    tuple val(sub), path("${sub}_orientation.npy"), emit: coords
    tuple val(sub), path("${sub}_history.txt"), emit: history

    shell:
    '''
    /scripts/optimize_fem.py !{msh} !{weights} !{centroid} \
                             !{coil} \
                             !{sub}_orientation.npy \
                             --out_msh !{sub}_optimized_fields.msh \
                             --out_geo !{sub}_optimized_coil.geo \
                             --history !{sub}_history.txt \
                             --ncores !{params.bayes_cpus.intdiv(2) - 2} \
                             grid \
                             --n_locations !{params.positional_grid_num} \
                             --n_rotations !{params.rotational_grid_num}
    '''
}

workflow optimize_wf{

    take:
       msh
       weights
       centroid
       coil

    main:
        i_grid_optimization = msh.join(weights)
                                .join(centroid)
                                .spread([coil])
        grid_optimization(i_grid_optimization)

    emit:
        history = grid_optimization.out.history
        coil = grid_optimization.out.coil
        fields = grid_optimization.out.fields
        coords = grid_optimization.out.coords
}

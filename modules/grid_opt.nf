nextflow.preview.dsl = 2

// Set up grid optimization schema
process grid_optimization{

    stageInMode 'copy'
    label 'rtms'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh),\
    path(weights), path(centroid), path(coil)

    output:
    tuple val(sub), path("${sub}_orientation.txt"), emit: orientation
    tuple val(sub), path("${sub}_history.txt"), emit: history

    shell:
    '''
    #!/bin/bash

    /scripts/gridopt_fem.py !{msh} !{weights} !{centroid} \
                            !{coil} \
                            $(pwd)/!{sub}_orientation.txt \
                            !{params.positional_grid_num} \
                            !{params.rotational_grid_num} \
                            --history !{sub}_history.txt \
                            --workdir $(pwd) \
                            --ncpus !{params.grid_cpus} \
                            --batchsize !{params.batch_size}
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
                                .spread([coil]) | view
        grid_optimization(i_grid_optimization)

    emit:
        orientation = grid_optimization.out.orientation
        history = grid_optimization.out.history
}

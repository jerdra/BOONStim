nextflow.preview.dsl = 2

// Set up optunity optimization schema
process optunity_optimization{

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

    /scripts/optunityopt_fem.py !{msh} !{weights} !{centroid} \
                            !{coil} \
                            $(pwd)/!{sub}_orientation.txt \
                            !{params.position_optunity_num} \
                            !{params.rotational_optunity_num} \
                            --history !{sub}_history.txt \
                            --workdir !(pwd)
                            --ncpus !{params.optunity_cpus} \
                            --batchsize !{params.batch_size}
    '''

}

workflow optimize{

    take:
       msh,
       weights,
       centroid,
       coil

    main:
        i_optunity_optimization = msh.join(weights)
                                .join(centroid)
                                .join(coil)
        optunity_optimization(i_optunity_optimization)

    emit:
        orientation = optimize_coil.out.orientation
        history = optimize_coil.out.history
}

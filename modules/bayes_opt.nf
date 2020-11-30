nextflow.preview.dsl = 2

// Process definitions
process bayesian_optimization{

    stageInMode 'copy'
    label 'rtms'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh), path(weights),\
          path(centroid), path(coil)

    output:
    tuple val(sub), path("${sub}_orientation.txt"), emit: orientation
    tuple val(sub), path("${sub}_history.txt"), emit: history

    shell:
    '''
    /scripts/optimize_fem.py !{msh} !{weights} !{centroid} \
                             !{coil} \
                             !{sub}_orientation.txt \
                             --history !{sub}_history.txt \
                             --n-iters !{params.max_iters} \
                             --cpus !{params.bayes_cpus.intdiv(2) - 2} \
                             --tmp-dir $(pwd)
    '''
}

workflow optimize_wf{

    take:
        msh
        weights
        centroid
        coil

    main:
        i_bayesian_optimization = msh
                                    .join(weights)
                                    .join(centroid)
                                    .join(coil)
        bayesian_optimization(i_bayesian_optimization)

    emit:
        orientation = bayesian_optimization.out.orientation
        history = bayesian_optimization.out.history
}

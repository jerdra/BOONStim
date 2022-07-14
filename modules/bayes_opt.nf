nextflow.preview.dsl = 2

// Process definitions
process bayesian_optimization{

    stageInMode 'copy'
    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh), path(weights),\
          path(centroid), path(coil)

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
                             --nworkers 4 \
                             --ncores !{params.bayes_cpus} \
                             bayesian \
                             --max_iterations !{params.max_iters}
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
                                    .spread([coil])
        bayesian_optimization(i_bayesian_optimization)

    emit:
        history = bayesian_optimization.out.history
        coil = bayesian_optimization.out.coil
        fields = bayesian_optimization.out.fields
        coords = bayesian_optimization.out.coords
}

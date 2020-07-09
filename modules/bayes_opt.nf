nextflow.preview.dsl = 2

// Process definitions
process optimize_coil{

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
                             --n-iters 50 \
                             --cpus 5 \
                             --tmp-dir /tmp
    '''
}

workflow optimize_wf{

    take:
        msh
        weights
        centroid
        coil

    main:
        i_optimize_coil = msh
                            .join(weights)
                            .join(centroid)
                            .join(coil)
        optimize_coil(i_optimize_coil)

    emit:
        orientation = optimize_coil.out.orientation
        history = optimize_coil.out.history
}

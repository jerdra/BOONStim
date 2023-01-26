nextflow.preview.dsl = 2

process bayesian_optimization{

    /*
    Bayesian optimization of E-field magnitude on a continuous-valued
    scalar map and realistic head model

    Arguments:
        sub (str): Subject ID
        msh (Path): Path to .msh file
        weights (Path): .npy scalar map
        centroid (Path): Path to seed coordinate
        coil (Path): Path to .nii.gz or .ccd coil file

    Parameters:
        bin (Path): Path to BOONStim scripts directory
        bayes_cpus (int): Number of threads to use for simulations
        max_iters (int): Maximum number of iterations allowed

    Outputs:
        fields (channel): (sub, msh: Path) Optimal E-field simulation .msh
        coil (channel): (sub, geo: Path) Optimal E-field coil position .geo
        coords (channel): (sub, coords: Path) Optimal coil orientation matrix
        history (channel): (sub, history: Path) Record of scores across iterations
    */


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

    /*
    Perform Bayesian Optimization on continuous-valued target map

    Arguments:
        msh (channel): (subject, msh) Subject msh files
        weights (channel): (subject, wf) Subject target scalar maps
        centroid (channel): (subject, centroid) Subject seed coordinates
        coil (value): .nii.gz or .ccd Coil dA/dt or definition file respectively

    Parameters:
        bin (Path): Path to BOONStim scripts directory
        bayes_cpus (int): Number of threads to use for simulations
        max_iters (int): Maximum number of iterations allowed

    Outputs:
        fields (channel): (sub, msh: Path) Optimal E-field simulation .msh
        coil (channel): (sub, geo: Path) Optimal E-field coil position .geo
        coords (channel): (sub, coords: Path) Optimal coil orientation matrix
        history (channel): (sub, history: Path) Record of scores across iterations
    */

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

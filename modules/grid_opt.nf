nextflow.preview.dsl = 2

// Set up grid optimization schema
process grid_optimization{

    /*
    Grid optimization of E-field magnitude on a continuous-valued
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
        positional_grid_num (int): Number of locations to sample 
            on each axis (a total of N X N locations will be sampled)
        rotational_grid_num (int): Number of orientations to sample per location

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
        positional_grid_num (int): Number of locations to sample 
            on each axis (a total of N X N locations will be sampled)
        rotational_grid_num (int): Number of orientations to sample per location

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

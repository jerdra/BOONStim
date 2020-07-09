nextflow.preview.dsl = 2


include optimize_wf as optimize from "${params.optimization_module}" params(params)

process evaluate_fem{

    input:
    tuple val(sub), path(msh),\
    path(weights), path(centroid),\
    path(orientation), path(coil)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"),\
    path("${sub}_optimized_coil.geo"),\
    path("${sub}_optimized_coords.npy"), emit: sim_results

    shell:
    '''
    import numpy as np
    from fieldopt.objective import FieldFunc

    coords = np.genfromtxt("!{orientation}")
    centroid = np.genfromtxt("!{centroid}")
    wf = np.load("!{weights}")

    fem = FieldFunc("!{msh}",
                    initial_centroid=centroid,
                    tet_weights=wf,
                    coil="!{coil}",
                    field_dir=!(pwd),
                    cpus=2
                   )

    _, matsimnibs = fem.run_simulation(coords, "!{sub}_optimized_fields.msh",
                        "!{sub}_optimized_coil.geo")
    np.save("!{sub}_optimized_coords.npy", matsimnibs)

    '''

}

workflow optimize_wf{

    take:
       msh
       weights
       centroid
       coil

    main:

        // Run optimization routine
        optimize(msh, weights, centroid, coil)

        // Evaluate optimal mesh
        i_evaluate_fem = msh.join(weights)
                            .join(centroid)
                            .join(optimize.out.orientation)
                            .spread([coil])
        evaluate_fem(i_evaluate_fem)

    emit:
        orientation = optimize.out.orientation
        history = optimize.out.history
        sim_results = evaluate_fem.out.sim_results
}

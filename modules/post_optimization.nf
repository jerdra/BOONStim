nextflow.preview.dsl = 2

process evaluate_fem{

    input:
    tuple val(sub), path(mesh), path(orientation), path(coil)

    output:
    tuple val(sub), path(sim_msh), path(coil_geo), emit: sim_result

    shell:
    '''
    #!/usr/bin/env python

    import numpy as np
    from fieldopt.objective import FieldFunc

    matsimnibs = np.load("!{orientation}")
    f = FieldFunc("!{mesh}", 
    '''
}

workflow post_optimization_wf{

    take:
        mesh
        orientation
        history
        coil


    main:

        // Run the optimal simulation
        i_evaluate_fem = mesh.join(orientation)
                             .spread([coil])
        evaluate_fem(i_evaluate_fem)

    emit:

}

nextflow.preview.dsl = 2


include optimize_wf as optimize from "${params.optimization_module}" params(params)

process evaluate_fem{

    label 'rtms'

    input:
    tuple val(sub), path(msh),\
    path(weights), path(centroid),\
    path(orientation), path(coil)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: fields
    tuple val(sub), path("${sub}_optimized_coil.geo"), emit: coil
    tuple val(sub), path("${sub}_optimized_coords.npy"), emit: coords

    shell:
    '''
    #!/usr/bin/env python
    import os
    import numpy as np
    from fieldopt.objective import FieldFunc

    coords = np.genfromtxt("!{orientation}")
    centroid = np.genfromtxt("!{centroid}")
    wf = np.load("!{weights}")

    fem = FieldFunc("!{msh}",
                    initial_centroid=centroid,
                    tet_weights=wf,
                    coil="!{coil}",
                    field_dir=os.getcwd(),
                    cpus=2
                   )

    _, matsimnibs = fem.run_simulation(coords, "!{sub}_optimized_fields.msh",
                        "!{sub}_optimized_coil.geo")

    # Save matsimnibs matrix
    np.save("!{sub}_optimized_coords.npy", matsimnibs)
    '''

}


process brainsight_transform{

    label 'rtms'

    input:
    tuple val(sub), path(orientation)

    output:
    tuple val(sub), path("${sub}_brainsight.csv"), emit: brainsight_coords

    shell:
    '''
    #!/usr/bin/env python

    from numpy import arcsin, degrees, arctan2
    import numpy as np

    matsimnibs = np.load("!{orientation}")

    matsimnibs[:3,2] = -msn[:3,2]
    matsimnibs[:3,0] = -msn[:3,0]

    xyz_alpha = degrees(arctan2(-msn[1,2],msn[2,2]))
    xyz_beta = degrees(arcsin(msn[0,2]))
    xyz_gamma = degrees(arctan2(-msn[0,1],msn[0,0,]))

    to_write = np.array([
        matsimnibs[0,3],
        matsimnibs[1,3],
        matsimnibs[2,3],
        -xyz_alpha,
        -xyz_beta,
        xyz_gamma
    ])

    header = "X,Y,Z,AP,LR,Twist"

    np.savetxt("!{sub}_brainsight.csv", to_write,
        delimiter=",", header=header)
    '''

}

process localite_transform{

    label 'rtms'

    input:
    tuple val(sub), path(orientation)

    output:
    tuple val(sub), path("${sub}_localite.csv"), emit: localite_coords

    shell:
    '''
    #!/usr/bin/env python

    import numpy as np

    AFFINE = np.array([
    [0, 0, 1, 0],
    [0, -1, 0, 0,],
    [1, 0, 0, 0],
    [0, 0, 0, 1]
    ])

    matsimnibs = np.load("!{orientation}")
    localite_af = AFFINE @ matsimnibs
    np.savetxt("!{sub}_localite.csv", localite_af, delimiter=',')
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
        brainsight_transform(evaluate_fem.out.coords)
        localite_transform(evaluate_fem.out.coords)

    emit:
        history = optimize.out.history
        fields = evaluate_fem.out.fields
        coil = evaluate_fem.out.coil
        localite = localite_transform.out.localite_coords
        brainsight = brainsight_transform.out.brainsight_coords
        matsimnibs = evaluate_fem.out.coords
}

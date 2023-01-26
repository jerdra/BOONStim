nextflow.preview.dsl = 2


include {optimize_wf as optimize} from "${params.optimization_module}" params(params)

process publish_opt{

    publishDir path: "${params.out}/boonstim/${sub}/optimization", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(fields), path(coil), path(history)

    output:
    tuple path(fields), path(coil), path(history)

    shell:
    '''
    #!/bin/bash
    echo "Transferring optimization results to boonstim/!{sub}/optimization..."
    '''
}

process qc_parameteric_surf{
    /*
    Generate QC image of the parameterized sampling domain
    over the head

    Arguments:
        sub (str): Subject key
        msh (Path): Path to .msh file
        centroid (Path): Path to seed coordinate
        dscalar (Path): scalar-values target map

    Outputs:
        qc_parameteric (channel): (sub, qc: Path) HTML QC file of parameteric surface
    */


    label 'fieldopt'

    input:
    tuple val(sub), path(msh),\
    path(centroid), path(dscalar)

    output:
    tuple val(sub), path("${sub}_parameteric_qc.html"), emit:qc_parameteric

    shell:
    '''
    #!/bin/bash
    /scripts/qc/parameteric_surf.py !{msh} !{centroid} !{dscalar} \
                                    !{sub}_parameteric_qc.html

    '''
}


process brainsight_transform{

    /*
    Compute BrainSight X, Y, Z, AP, LR, Twist coordinates
    from a coil orientation matrix

    Arguments:
        sub (str): Subject ID
        orientation (Path): Coil orientation matsimnibs matrix

    Outputs:
        brainsight_coords (channel): (sub, brainsight: Path) Brainsight coordinates
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(orientation)

    output:
    tuple val(sub), path("${sub}_brainsight.csv"), emit: brainsight_coords

    shell:
    '''
    #!/usr/bin/env python

    from numpy import arcsin, degrees, arctan2
    import numpy as np

    msn = np.load("!{orientation}")

    msn[:3,1] = -msn[:3,1]
    msn[:3,2] = -msn[:3,2]

    xyz_alpha = degrees(arctan2(-msn[1,2],msn[2,2]))
    xyz_beta = degrees(arcsin(msn[0,2]))
    xyz_gamma = degrees(arctan2(-msn[0,1],msn[0,0]))

    to_write = np.array([
        msn[0,3],
        msn[1,3],
        msn[2,3],
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
    /*
    Compute Localite affine matrix from a coil orientation matrix

    Arguments:
        sub (str): Subject ID
        orientation (Path): Coil orientation matsimnibs matrix

    Outputs:
        localite_coords (channel): (sub, localite: Path) Localite affine matrix
    */

    label 'fieldopt'

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
    localite_af = matsimnibs @ AFFINE
    np.savetxt("!{sub}_localite.csv", localite_af, delimiter=',')
    '''
}

workflow optimize_wf{

    /*
    Perform TMS coil position optimization over a scalar-valued target map

    Arguments:
        msh (channel): (subject, msh: Path) .msh file
        weights (channel): (subject, wf) Subject target scalar maps
        centroid (channel): (subject, centroid) Subject seed coordinates
        coil (value): .nii.gz or .ccd Coil dA/dt or definition file respectively

    Outputs:
        fields (channel): (sub, msh: Path) Optimal E-field simulation .msh
        coil (channel): (sub, geo: Path) Optimal E-field coil position .geo
        coords (channel): (sub, coords: Path) Optimal coil orientation matrix
        history (channel): (sub, history: Path) Record of scores across iterations
        localite (channel): (sub, localite: Path) Localite affine matrix
        brainsight (channel): (sub, brainsight: Path) BrainSight coordinates
        matsimnibs (channel): (sub, matsimnibs: Path) Optimal matsimnibs coil orientation matrix
        qc_parsurf (channel): (sub, qchtml: Path) Parameteric domain QC page
    */


    take:
       msh
       weights
       centroid
       coil

    main:

        // Run optimization routine
        optimize(msh, weights, centroid, coil)

        brainsight_transform(optimize.out.coords)
        localite_transform(optimize.out.coords)
        i_qc_parameteric_surf = msh.join(centroid)
                                   .join(weights)
        qc_parameteric_surf(i_qc_parameteric_surf)

        i_publish_out = optimize.out.fields
            .join(optimize.out.coil)
            .join(optimize.out.history)
        publish_opt(i_publish_out)


    emit:
        history = optimize.out.history
        fields = optimize.out.fields
        coil = optimize.out.coil
        matsimnibs = optimize.out.coords
        qc_parsurf = qc_parameteric_surf.out.qc_parameteric
}

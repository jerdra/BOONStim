nextflow.preview.dsl=2

process reorient_coil {

    label 'numpy'

    /*
    Re-orient coil so that handle is pointing posterior

    Arguments:
        sub (str): Subject ID
        msn (Path): Matsimnibs matrix

    Outputs:
        reoriented (channel): (sub, reoriented: Path) Matsimnibs reoriented so that
            coil handle is pointing posterior
        flipped (channel): (sub, flipped: Path)  Whether the coil handle was flipped 
            or not
    */
    input:
    tuple val(sub), path(msn)

    output:
    tuple val(sub), path("${sub}_reoriented.npy"), emit: reoriented
    tuple val(sub), path("${sub)_isflip.txt"), emit: flipped

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    msn = np.load('${msn}')
    coil_anterior_vec = msn[:3, 1]
    isflip = 0

    # anterior component of coil head points posterior
    if msn[1, 1] < 0:
        msn[:, 1] *= -1

        # Maintain right-handedness
        msn[:, 0] *= -1
        isflip = 1

    msn.save('${sub}_reoriented.npy')
    with open('${sub}_isflip.txt', 'w') as f:
        f.write(str(isflip))
    """
        
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

workflow neuronav_wf {

    /*
    Produce Neuronavigation (BrainSight, Localite) outputs
    from SimNIBS coil orientation matrix

    Arguments:
        msn (channel): (subject, msn: Path) Optimal matsimnibs matrix

    Outputs:
        brainsight (channel): (subject, brainsight: Path) Brainsight coordinates
        localite (channel): (subject, localite: Path) Localite matrix
        isflipped (channel): (subject, flip: bool) Should current be flipped?

    Note:
        `flip` will be set to true if a coil orientation flip (180 deg) is required
        to ensure that the coil handle is pointing in the posterior direction
    */

    take:
        msn

    main:
        reorient_coil(msn)
        brainsight_transform(reorient_coil.out.reoriented)
        localite_transform(reorient_coil.out.reoriented)

    emit:
        brainsight = brainsight_transform.out.brainsight_coords
        localite = localite_transform.out.localite_coords
        isflipped = reorient_coil.out.flipped
}

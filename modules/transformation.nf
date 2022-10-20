nextflow.preview.dsl=2

process rigidRegistration {

    /*
    Perform rigid registration between two images

    Arguments:
        sub (str): Subject ID
        moving (Path): Path to source image
        fixed (Path): Path to reference image

    Outputs:
        transformed (channel): (sub, transformed) Moving image transformed to Fixed
        transform (channel): (sub, transform) Rigid-body transformation from moving
            to fixed
    */

    label 'ants'

    input:
    tuple val(sub), path(moving), path(fixed)

    output:
    tuple val(sub), path("${sub}_transformed.nii.gz"), emit: transformed
    tuple val(sub), path("${sub}_transformation.mat"), emit: transform

    shell:
    """
    antsRegistration --dimensionality 3 --float 0 \
    --interpolation Linear \
    --winsorize-image-intensities [0.005, 0.995] \
    --use-histogram-matching 0 \
    --transform Rigid[0.1] \
    --metric MI[$moving,$fixed,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --output [${sub}_transformed_,${sub}_transformed.nii.gz]

    mv ${sub}_transformed_0GenericAffine.mat ${sub}_transformation.mat

    """
}

process coordinate_transform {

    /*
    Apply ANTS coordinate transform on input
    coordinates
    Ensure that transform being used is the forward image transformed
    The inverse flag is set in this function

    Arguments:
        sub (str): Subject ID
        coords (Path): input coordinates
        transform (Path): Ants transformation file

    Outputs:
        transformed (channel): (sub, transformed: Path)
    */

    label 'ants'

    input:
    tuple val(sub), path(coords), path(transform)

    output:
    tuple val(sub), path("${sub}_target_coord.csv"), emit: transformed

    shell:
    """
    antsApplyTransformsToPoints \
        -d 3 \
        -i ${coords} \
        -t [${transform},0] \
        -o ${sub}_target_coord.csv
    """
}


process format_for_ants {

    /*
    Format a RAS coordinate into an LPS ANTS CSV file

    Arguments:
        sub (str): Subject ID
        coords (Path): Path to coordinate .txt file

    Output:
        ants_csv (channel): (sub, csv: Path) Ants formatted CSV file
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_pretransform.csv"), emit: ants_csv

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    # Read in RAS, convert to LPS
    coords = np.loadtxt('${coords}', delimiter=' ')
    coords[0] = -coords[0]
    coords[1] = -coords[1]

    # Save with ANTS header
    header = 'x,y,z,t,label'
    coords_as_txt = ','.join([str(x) for x in coords])

    # ANTS wants a label column
    coords_as_txt += ",0,nolabel"

    with open('${sub}_pretransform.csv', 'w') as f:
        f.writelines([header, '\\n', coords_as_txt, '\\n'])
    """
}

process format_from_ants {

    /*
    Format an  LPS ANTS CSV file into a RAS coordinate file

    Arguments:
        sub (str): Subject ID
        coords (Path): Path to ants coordinate .txt file

    Output:
        ras_coords (channel): (sub, npy: Path) Path to .npy file in RAS
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_ras_coords.npy"), emit: ras_coords

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    # LPS space
    coords = np.loadtxt('${coords}', skiprows=1, usecols=(0,1,2,3), delimiter=',')
    coords[0] = -coords[0]
    coords[1] = -coords[1]
    np.save("${sub}_ras_coords.npy", coords[:-1])
    """
}

workflow map_coordinate {
    /*
    Perform rigid-body coordinate mapping from a set of RAS coordinates
    using an ANTS transformation

    Arguments:
        coordinates (channel): (sub, coordinate: Path) RAS coordinate files
        transforms (channel): (sub, transform: Path) ANTS Rigid-body transformtion

    Outputs:
        transformed (channel): (sub, transformed: Path) RAS coordinates in target space

    */

    take:
        coordinates
        transforms

    main:
        format_for_ants(coordinates)
        coordinate_transform(
            format_for_ants.out.ants_csv.join(transforms)
        )
        format_from_ants(coordinate_transform.out.transformed)

    emit:
        transformed = format_from_ants.out.ras_coords

}

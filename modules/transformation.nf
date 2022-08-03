include { rigidRegistration, coordinate_transform } from "../modules/utils.nf"

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
    --output [${sub}_transformed_,warped.nii.gz]

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

    input:
    tuple val(sub), path(coords), path(transform)

    output:
    tuple val(sub), path("${sub}_target_coord.txt"), emit: transformed

    shell:
    """
    antsApplyTransformToPoints \
        -d 3 -e 0 \
        -i ${coords} \
        -t [${transform},1] \
        -o ${sub}_target_coord.txt
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

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_pretransform.txt"), emit: ants_csv

    shell:
    """
    echo "x,y,z,t" > "${sub}_pretransform.txt"
    echo "${x},${y},${z},0" >> "${sub}_pretransform.txt"
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

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_ras_coords.npy"), emit: ras_coords

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    # LPS space
    coords = np.loadtxt('${sub}_ras_coords.txt')
    coords[0] = -coords[0]
    coords[1] = -coords[1]
    np.savetxt("${sub}_ras_coords.txt", coords[:-1], delim=',')
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
        coordinate_transform(coordinates.join(transforms))
        format_from_ants(coordinate_transform.out.transformed)

    emit:
        transformed = format_from_ants.ras_coords

}

nextflow.preview.dsl=2

process apply_mask {

    /*
    Apply a mask to a dscalar image

    Arguments:
        sub (str): Subject ID
        dscalar (Path): Input dscalar file to mask
        mask (Path): Mask dscalar file

    Outputs:
        masked (channel): (sub, masked: Path)
    */

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(mask)

    output:
    tuple val(sub), path("${sub}.masked.dscalar.nii"), emit: masked

    shell:
    '''
    wb_command -cifti-math \
                "x * (mask > 0)" \
                -var "x" !{dscalar} \
                -var "mask" !{mask} \
                !{sub}.masked.dscalar.nii
    '''
}

process cifti_dilate {

    /*
    Dilate an input dscalar file

    Arguments:
        sub (str): Subject ID
        dscalar (Path): Path to dscalar file
        left (Path): Path to left surface file (midthickness typically)
        right (Path): Path to right surface file (midthickness typically)

    Outputs:
        dilated (channel): (sub, dilated: Path)
    */

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(left), path(right)

    output:
    tuple val(sub), path("${sub}.dilated.dscalar.nii"), emit: dilated

    shell:
    '''
    wb_command -cifti-dilate \
                !{dscalar} \
                COLUMN \
                6 6 \
                -left-surface !{left} \
                -right-surface !{right} \
                !{sub}.dilated.dscalar.nii
    '''

}

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

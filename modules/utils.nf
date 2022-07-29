nextflow.preview.dsl=2

process apply_mask {

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
    --metric MI[$t1brain,$template,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --output [${sub}_transformed_,warped.nii.gz]

    mv ${sub}_transformed_0GenericAffine.mat ${sub}_transformation.mat

    """
}

process coordinate_transform {

    /*
    * Apply ANTS coordinate transform on input
    * coordinates
    * Ensure that transform being used is the forward image transformed
    * The inverse flag is set in this function
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

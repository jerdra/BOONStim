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

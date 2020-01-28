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
                -var "x" !{weightfunc} \
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
                ROW \
                6 6 \
                !{sub}.dilated.scalar.nii
    '''

}

nextflow.preview.dsl=2


process project_mask2surf{

    label 'connectome'
    input:
    tuple val(sub), path(white), path(pial), path(midthick), path(mask)

    output:
    tuple val(sub), path('surfmask.shape.gii'), emit: mask_shape

    shell:
    '''
    wb_command -volume-to-surface-mapping \
                !{mask} \
                !{midthick} \
                -ribbon-constrained \
                    !{white} !{pial} \
                surfmask.shape.gii
    '''

}

process binarize_mask{

    label 'connectome'
    input:
    tuple val(sub), path(surfmask)
    
    output:
    tuple val(sub), path('bin_mask.shape.gii'), emit:bin_mask

    shell:
    '''
    wb_command -metric-math \
                "(x>0)" \
                -var "x" !{surfmask} \
                bin_mask.shape.gii
    '''

}

process dilate_mask{

    label 'connectome'
    input:
    tuple val(sub), path(surfmask), path(midthick)

    output:
    tuple val(sub), path('dilated_mask.shape.gii'), emit:dilated_mask

    shell:
    '''
    wb_command -metric-dilate \
                !{surfmask} \
                !{midthick} \
                12 -nearest \
                dilated_mask.shape.gii
    '''
}

process split_weightfunc {

    label 'connectome'
    input:
    tuple val(sub), path(weightfunc)

    output:
    tuple val(sub), path('weightfunc.L.shape.gii'), emit:left_weightfunc
    tuple val(sub), path('weightfunc.R.shape.gii'), emit:right_weightfunc

    shell:
    '''
    wb_command -cifti-separate \
                !{weightfunc} \
                COLUMN \
                -metric CORTEX_LEFT weightfunc.L.shape.gii \
                -metric CORTEX_RIGHT weightfunc.R.shape.gii
                
    '''

}

process apply_mask {

    label 'connectome'
    input:
    tuple val(sub), path(weightfunc), path(mask)
    
    output:
    tuple val(sub), path('masked_weightfunc.L.shape.gii'), emit:masked_weightfunc_shape

    shell:
    '''
    wb_command -metric-mask \
                !{weightfunc} \
                !{mask} \
                masked_weightfunc.L.shape.gii
    '''
}

process threshold_weightfunc {

    label 'connectome'
    input:
    tuple val(sub), path(weightfunc)
    
    output:
    tuple val(sub), path('thresholded_masked_weightfunc.L.shape.gii'), emit:thresholded_weightfunc_shape

    shell:
    '''
    wb_command -metric-math \
            "x * (x>0)" \
            -var x !{weightfunc} \
            thresholded_masked_weightfunc.L.shape.gii
    '''
}

process recombine_weightfunc{

    label 'connectome'
    input:
    tuple val(sub), path(left_shape), path(right_shape)

    output:
    tuple val(sub), path('masked_weightfunc.dscalar.nii'), emit: weighted_mask
    
    shell:
    '''
    wb_command -cifti-create-dense-scalar \
                masked_weightfunc.dscalar.nii \
                -left-metric !{left_shape} \
                -right-metric !{right_shape}
    '''

}

//RECOMBINE INTO DSCALAR

workflow mask_wf {

    get: 
        cifti
        weightfile

    main:
        project_mask_inputs = cifti
                                .map{ s,c ->    [
                                                    s,
                                                    "${c}/${s}/MNINonLinear/fsaverage_LR32k/${s}.L.white.32k_fs_LR.surf.gii",
                                                    "${c}/${s}/MNINonLinear/fsaverage_LR32k/${s}.L.pial.32k_fs_LR.surf.gii",
                                                            "${c}/${s}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                                    "${params.mask}"
                                                ]
                                    }

        // Project, binarize, dilate
        project_mask2surf(project_mask_inputs)
        binarize_mask(project_mask2surf.out.mask_shape)

        // Add midthickness file
        dilate_mask_input = binarize_mask.out.bin_mask
                                        .join(cifti, by: 0)
                                        .map{ s,b,c ->  [   
                                                            s,b,
                                                            "${c}/${s}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii"
                                                        ]
                                            }
        dilate_mask(dilate_mask_input)
        dilate_mask.out.dilated_mask

        // Split weightfunction
        split_weightfunc(weightfile)

        // Pair up weightfunc --> mask and apply
        apply_mask_inputs = split_weightfunc.out.left_weightfunc
                                            .join(dilate_mask.out.dilated_mask, by: 0)
        apply_mask(apply_mask_inputs)

        // Threshold
        threshold_weightfunc(apply_mask.out.masked_weightfunc_shape)

        // Recombine weightfunction
        recombine_inputs = threshold_weightfunc.out.thresholded_weightfunc_shape
                                            .join(split_weightfunc.out.right_weightfunc, by: 0)
        recombine_weightfunc(recombine_inputs)

    emit:
        weighted_mask = recombine_weightfunc.out.weighted_mask

}

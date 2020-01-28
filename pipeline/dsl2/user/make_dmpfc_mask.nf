nextflow.preview.dsl=2


process project_mask2surf{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(white), path(pial), path(midthick), path(mask)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.roi_groupmask.shape.gii"), emit: mask_shape

    shell:
    '''
    wb_command -volume-to-surface-mapping \
                !{mask} \
                !{midthick} \
                -ribbon-constrained \
                    !{white} !{pial} \
                !{sub}.!{hemi}.roi_groupmask.shape.gii
    '''

}

process binarize_mask{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(surfmask)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.bin_mask.shape.gii"), emit:bin_mask

    shell:
    '''
    wb_command -metric-math \
                "(x>0)" \
                -var "x" !{surfmask} \
                !{sub}.!{hemi}.bin_mask.shape.gii
    '''

}

process dilate_mask{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(surfmask), path(midthick)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.dilated_mask.shape.gii"), emit:dilated_mask

    shell:
    '''
    wb_command -metric-dilate \
                !{surfmask} \
                !{midthick} \
                12 -nearest \
                !{sub}.!{hemi}.dilated_mask.shape.gii
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
    tuple val(sub), path("${sub}.masked_weightfunc.dscalar.nii"), emit:masked_weightfunc

    shell:
    '''
    wb_command -cifti-math \
                "x * (mask > 0)" \
                -var "x" !{weightfunc} \
                -var "mask" !{mask} \
                !{sub}.masked_weightfunc.dscalar.nii
    '''
}

process threshold_weightfunc {

    label 'connectome'
    input:
    tuple val(sub), path(weightfunc)

    output:
    tuple val(sub), path("${sub}.thresholded_masked_weightfunc.dscalar.nii"), emit:thresholded_weightfunc

    shell:
    '''
    wb_command -cifti-math \
            "x * (x>0)" \
            -var x !{weightfunc} \
            !{sub}.thresholded_masked_weightfunc.dscalar.nii
    '''
}

process recombine_masks{

    label 'connectome'
    input:
    tuple val(sub), path(left), path(right)

    output:
    tuple val(sub), path("${sub}.weightfunc_mask.dscalar.nii"), emit: dscalar_mask

    shell:
    '''
    wb_command -cifti-create-dense-scalar \
                !{sub}.weightfunc_mask.dscalar.nii \
                -left-metric !{left} \
                -right-metric !{right}
    '''


}

//RECOMBINE INTO DSCALAR

workflow mask_wf {

    get:
        cifti
        weightfile

    main:
        project_lmask_inputs = cifti
                                .map{ s,c ->    [
                                                    s,
                                                    'L',
                                                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.white.32k_fs_LR.surf.gii",
                                                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.pial.32k_fs_LR.surf.gii",
                                                            "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                                    "${params.mask}"
                                                ]
                                    }

        project_rmask_inputs = cifti
                                .map{ s,c ->    [
                                                    s,
                                                    'R',
                                                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.white.32k_fs_LR.surf.gii",
                                                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.pial.32k_fs_LR.surf.gii",
                                                            "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.midthickness.32k_fs_LR.surf.gii",
                                                    "${params.mask}"
                                                ]
                                    }

        // Project mask to surfaces
        project_mask_inputs = project_lmask_inputs.mix(project_rmask_inputs)
        project_mask2surf(project_mask_inputs)

        // Binarize masks on surface
        binarize_mask(project_mask2surf.out.mask_shape)
        binarize_mask.out.bin_mask

        // Set up dilation input
        dilate_mask_input = binarize_mask.out.bin_mask
                                        .combine(cifti, by: 0)
                                        .map{ s,h,b,c ->  [
                                                            s,h,b,
                                                            "${c}/MNINonLinear/fsaverage_LR32k/${s}.${h}.midthickness.32k_fs_LR.surf.gii"
                                                        ]
                                            }
        dilate_mask(dilate_mask_input)

        // Recombine shape mask files
        recombine_input = dilate_mask.out.dilated_mask
                                        .map { s,h,d -> [s,d] }
                                        .groupTuple ( by: 0, sort: {it.baseName}, size: 2)
                                        .map { s,f -> [s, f[0], f[1]] }

        recombine_masks(recombine_input)
        recombine_masks.out.dscalar_mask

        // Apply bilateral mask
        apply_mask_input = weightfile
                                .join(recombine_masks.out.dscalar_mask, by:0)
        apply_mask(apply_mask_input)

        // Threshold
        threshold_weightfunc(apply_mask.out.masked_weightfunc)

    emit:
        weighted_mask = threshold_weightfunc.out.thresholded_weightfunc

}

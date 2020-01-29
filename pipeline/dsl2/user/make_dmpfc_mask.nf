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

process dilate_dmpfc_mask{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(dmpfc_mask), path(midthick)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.dmpfc_mask_dilated.shape.gii"), emit: dmpfc_mask

    shell:
    '''
    wb_command -metric-dilate \
                !{dmpfc_mask} \
                !{midthick} \
                6 \
                !{sub}.!{hemi}.dmpfc_mask_dilated.shape.gii
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

process make_symmetric_dscalar{

    label 'connectome'
    input:
    tuple val(sub), path(left_mask), path(right_mask)

    output:
    tuple val(sub), path("${sub}.dmpfc_mask_symmetric.dscalar.nii"), emit: dmpfc_mask

    shell:
    '''

    wb_command -metric-math \
                "(x+y) > 0" \
                -var x !{left_mask} \
                -var y !{right_mask} \
                !{sub}.L.dmpfc_symmetric.shape.gii

    wb_command -metric-math \
                "(x+y) > 0" \
                -var x !{right_mask} \
                -var y !{left_mask} \
                !{sub}.R.dmpfc_symmetric.shape.gii

    wb_command -cifti-create-dense-scalar \
                !{sub}.dmpfc_mask_symmetric.dscalar.nii \
                -left-metric !{sub}.L.dmpfc_symmetric.shape.gii \
                -right-metric !{sub}.R.dmpfc_symmetric.shape.gii
    '''
}


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

        // Project, binarize and dilate a priori mask
        project_mask_inputs = project_lmask_inputs.mix(project_rmask_inputs)
        project_mask2surf(project_mask_inputs)
        binarize_mask(project_mask2surf.out.mask_shape)
        dilate_dmpfc_mask(binarize_mask.out.bin_mask)

        // Make DMPFC symmetric mask
        make_symmetric_input = dilate_dmpfc_mask.out.dmpfc_mask
                                            .groupTuple(by: 0, size: 2, sort: {it})
                                            .map{ s,h,m -> [s, m[0], m[1]] }
        make_symmetric_dscalar(make_symmetric_input)

    emit:
        mask = make_symmetric_dscalar.out.dmpfc_mask

}

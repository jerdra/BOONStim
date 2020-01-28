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

process threshold_mask{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(surfmask)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.thres_mask.shape.gii"), emit: thres_mask

    shell:
    '''
    wb_command -metric-math \
                "(x>0.75)" \
                -var "x" !{surfmask} \
                !{sub}.!{hemi}.thres_mask.shape.gii
    '''
}

process pull_dmpfc_cluster{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(thres), path(midthick)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.dmpfc_mask.shape.gii"), emit: clusts

    shell:
    '''
    # Spatial clustering
    wb_command -metric-find-clusters \
                !{midthick} \
                !{thres} \
                0.5 50 \
                !{sub}.!{hemi}.thres_mask_clust.shape.gii

    # Pull maximal cluster # which will be dorsal DMPFC
    MAX_VAL=$(wb_command -metric-stats -reduce MAX !{sub}.!{hemi}.thres_mask_clust.shape.gii)

    wb_command -metric-math \
                "(x == round($MAX_VAL))" \
                -var "x" !{sub}.!{hemi}.thres_mask_clust.shape.gii \
                !{sub}.!{hemi}.dmpfc_mask.shape.gii
    '''
}

process dilate_cluster{

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

        // Threshold masks
        threshold_mask(project_mask2surf.out.mask_shape)

        // Pull DMPFC cluster
        dmpfc_cluster_input = threshold_mask.out.thres_mask
                                    .combine(cifti, by: 0)
                                    .map{ s,h,t,c ->[
                                                        s,h,t,
                                                        "${c}/MNINonLinear/fsaverage_LR32k/${s}.${h}.midthickness.32k_fs_LR.surf.gii"
                                                    ]
                                        }
        pull_dmpfc_cluster(dmpfc_cluster_input)

        // Dilate clusters
        dilate_cluster_input = pull_dmpfc_cluster.out.clusts
                                    .combine(cifti, by: 0)
                                    .map{ s,h,m,c ->    [
                                                            s,h,m,
                                                            "${c}/MNINonLinear/fsaverage_LR32k/${s}.${h}.midthickness.32k_fs_LR.surf.gii"
                                                        ]
                                        }
        dilate_cluster(dilate_cluster_input)

        // Make DMPFC mask symmetric
        make_symmetric_input = dilate_cluster.out.dmpfc_mask
                                            .groupTuple(by: 0, size: 2, sort: {it})
                                            .map{ s,h,m -> [s, m[0], m[1]] }
        make_symmetric_dscalar(make_symmetric_input)

    emit:
        mask = make_symmetric_dscalar.out.dmpfc_mask

}

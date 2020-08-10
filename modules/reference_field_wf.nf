nextflow.preview.dsl=2

process calc_distmap_from_coord{

    label 'connectome'

    input:
    tuple val(sub), val(hemi),\
    val(x), val(y), val(z),\
    path(surf)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.distmap.shape.gii"), emit: dist_shape

    shell:
    '''
    #Convert surfs --> metric
    wb_command -surface-coordinates-to-metric !{surf} surf.shape.gii

    # Subtract component
    wb_command -metric-math "x-!{x}" -var x surf.shape.gii -column "x coordinate" xdist.shape.gii
    wb_command -metric-math "x-!{y}" -var x surf.shape.gii -column "y coordinate" ydist.shape.gii
    wb_command -metric-math "x-!{z}" -var x surf.shape.gii -column "z coordinate" zdist.shape.gii

    # Compute euclidean distance
    wb_command -metric-math \
        "sqrt( x^2 + y^2 + z^2 )" \
        -var x xdist.shape.gii \
        -var y ydist.shape.gii \
        -var z zdist.shape.gii \
        !{sub}.!{hemi}.distmap.shape.gii
    '''
}

process join_distmaps{

    label 'connectome'

    input:
    tuple val(sub), path(L), path(R)

    output:
    tuple val(sub), path("${sub}.distmap.dscalar.nii"), emit: dist_dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar \
        !{sub}.distmap.dscalar.nii \
        -left-metric !{L} \
        -right-metric !{R}
    '''

}


process threshold_distmap {

    label 'connectome'

    input:
    tuple val(sub), path(distmap), val(thres)

    output:
    tuple val(sub), path("${sub}.distmap_roi.dscalar.nii"), emit: distmap_roi

    shell:
    '''
    wb_command -cifti-math \
        "x < !{thres}" \
        -var x !{distmap} \
        !{sub}.distmap_roi.dscalar.nii
    '''

}

process remove_medial_wall{

    label 'connectome'

    input:
    tuple val(sub), path(infile), path(tplfile)

    output:
    tuple val(sub), path("${sub}.medial_removed.dscalar.nii"), emit: cropped

    shell:
    '''
    #!/bin/bash
    wb_command -cifti-create-dense-from-template \
        !{tplfile} -cifti !{infile} \
        !{sub}.medial_removed.dscalar.nii
    '''
}

process get_precentral{

    label 'connectome'

    input:
    tuple val(sub), path(aparc)

    output:
    tuple val(sub), path("${sub}.precentral.dscalar.nii"), emit: precentral

    shell:
    '''
    #!/bin/bash

    wb_command -cifti-label-to-roi !{aparc} \
        -name L_precentral \
        left_precentral.dscalar.nii

    wb_command -cifti-label-to-roi !{aparc} \
        -name R_precentral \
        right_precentral.dscalar.nii

    wb_command -cifti-math "(x+y) > 0" \
        -var "x" left_precentral.dscalar.nii \
        -var "y" right_precentral.dscalar.nii \
        !{sub}.precentral.dscalar.nii
    '''
}

process dilate_mt_roi{

    label 'connectome'

    input:
    tuple val(sub), path(left_surf), path(right_surf), path(roi)

    output:
    tuple val(sub), path("${sub}.distmap_roi_dilated.dscalar.nii"), emit: distmap_roi

    shell:
    '''
    #!/bin/bash

    wb_command -cifti-dilate \
        !{roi} COLUMN 6 6 \
        -left-surface !{left_surf} \
        -right-surface !{right_surf} \
        -nearest \
        !{sub}.distmap_roi_dilated.dscalar.nii
    '''
}

process apply_precentral{

    label 'connectome'

    input:
    tuple val(sub), path(distmap), path(precentral)

    output:
    tuple val(sub), path("${sub}.distmap_precentral.dscalar.nii"), emit: distmap_roi

    shell:
    '''
    #!/bin/bash

    wb_command -cifti-math "x*y" \
        -var "x" !{distmap} \
        -var "y" !{precentral} \
        !{sub}.distmap_precentral.dscalar.nii
    '''
}


workflow calculate_reference_field_wf{

    take:
        ciftify

    main:

    i_calc_distmap_from_coord = ciftify
    .combine(["L","R"])
    .map{ s,c,h ->
        [
            s,h,
            "${params.ref_coord[0]}","${params.ref_coord[1]}","${params.ref_coord[2]}",
            "${c}/MNINonLinear/fsaverage_LR32k/${s}.${h}.pial.32k_fs_LR.surf.gii"
        ]
    }
    calc_distmap_from_coord(i_calc_distmap_from_coord)

    i_join_distmaps = calc_distmap_from_coord.out.dist_shape
                            .map{s,h,g -> [s,g]}
                            .groupTuple(by: 0, sort: { it.getBaseName() })
                            .map{s,g -> [s, g[0], g[1]]}
    join_distmaps(i_join_distmaps)

    i_threshold_distmap = join_distmaps.out.dist_dscalar
                            .map{s,d -> [s,d,"$params.ref_dist"]}
    threshold_distmap(i_threshold_distmap)

    // Get precentral sulcus
    i_get_precentral = ciftify
                        .map{ s,c ->
                        [
                            s,
                            "${c}/MNINonLinear/fsaverage_LR32k/${s}.aparc.32k_fs_LR.dlabel.nii"
                        ]}
    get_precentral(i_get_precentral)

    // Dilate the thresholded distance map to increase coverage
    i_dilate_mt = ciftify
                    .map{s,c ->
                    [
                    s,
                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                    "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.midthickness.32k_fs_LR.surf.gii"

                    ]}
                    .join(threshold_distmap.out.distmap_roi)

    dilate_mt_roi(i_dilate_mt)

    // Remove medial wall
    i_remove_medial_wall = dilate_mt_roi.out.distmap_roi
                                .join(i_get_precentral)
    remove_medial_wall(i_remove_medial_wall)

    // Apply precentral mask
    i_apply_precentral = remove_medial_wall.out.cropped
                            .join(get_precentral.out.precentral)
    apply_precentral(i_apply_precentral)



    emit:
        roi = apply_precentral.out.distmap_roi

}



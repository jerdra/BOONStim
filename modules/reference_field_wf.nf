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
            "${c}/MNINonLinear/fsaverage_LR32k/${s}.${h}.midthickness.32k_fs_LR.surf.gii"
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


    emit:
        roi = threshold_distmap.out.distmap_roi

}



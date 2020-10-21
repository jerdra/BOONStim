nextflow.preview.dsl = 2

process split_dscalar{

    label 'connectome'

    input:
    tuple val(id), path(dscalar)

    output:
    tuple val(id), val('L'), path('L.shape.gii'), emit: left
    tuple val(id), val('R'), path('R.shape.gii'), emit: right

    shell:
    '''
    wb_command -cifti-separate \
                !{dscalar} \
                COLUMN \
                -metric CORTEX_LEFT L.shape.gii \
                -metric CORTEX_RIGHT R.shape.gii
    '''

}

process resample_surf{

    label 'connectome'

    input:
    tuple val(id), val(hemi), path(shape), path(source_sphere), path(target_sphere)

    output:
    tuple val(id), val(hemi), path("resampled_shape.${hemi}.shape.gii"), emit: resampled

    shell:
    '''
    wb_command -metric-resample \
                !{shape} \
                !{target_sphere} \
                !{source_sphere} \
                BARYCENTRIC \
                resampled_shape.!{hemi}.shape.gii
    '''


}

process recombine {

    label 'connectome'

    input:
    tuple val(id), path(left), path(right)

    output:
    tuple val(id), path("resampled_dscalar.nii"), emit: dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar \
                 resampled_dscalar.nii \
                -left-metric !{left} \
                -right-metric !{right}
    '''


}

// Workflow to take fsaverage_LR32k surfaces and resample them to native space
workflow resample2native_wf {

    take:
        dscalar
        msm_sphere

    main:

        /* Associate each dscalar map with a unique ID, since we may
        have more than 1 map/subject */
        udscalar = dscalar.map{a ->
                               [ a[0], UUID.randomUUID().toString(),
                                 a[1..<-1], a[-1] ]}
                          .multiMap{s, u, ids, d ->
                            map2id: [u, s, ids]
                            map2dscalar: [u, d]
                            map2sub: [u, s]}

        split_dscalar(udscalar.map2dscalar)

        // Map back to subject, then append associated spheres
        resample_input=\
            split_dscalar.left.mix(split_dscalar.right)
                              .join(udscalar.map2sub, by: 0)
                              .map{u, h, d, s -> [s, h, u, d]}
                              .combine(msm_sphere, by: [0,1])
                              .map{s,h,u,d,m ->
                                    [
                                        u,h,d,m,
                                        "${params.atlas}/${h}.sphere.32k_fs_LR.surf.gii"
                                    ]}
        resample_surf(resample_input)

        // Recombine based on UUID assigned to original dscalar
        recombine_input = resample_surf.out.resampled
                                        .map { u,h,f -> [u,f] }
                                        .groupTuple ( by: 0, sort: { it.baseName }, size: 2 )
                                        .map { u,f -> [u,f[0],f[1]] }

        recombine(recombine_input)

        // Re-form the original structure of the input
        output = udscalar.map2id.join(recombine.out.dscalar, by: 0)
                        .map{u, s, ids, d -> [s, ids, d].flatten()}


        emit:
            resampled = output

}

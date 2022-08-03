nextflow.preview.dsl = 2

process split_dscalar{

    /*
    Split a surface-based dscalar file into .shape.gii files

    Arguments:
        id (str): ID key
        dscalar (Path): dscalar file

    Outputs:
        left (channel): (id, 'L', left: Path): Left GIFTI
        righL (channel): (id, 'R', left: Path): Right GIFTI

    Note: Subcortical volume will be discarded in this step
    */

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

    /*
    Perform surface-based resampling

    Arguments:
        id (str): ID Key
        hemi (Union['L', 'R']): Hemisphere key
        shape (Path): GIFTI shape file to resample
        source_sphere (Path): Sphere associated with GIFTI shape file (fsaverage_LR32k)
        target_sphere (Path): Sphere to resample to

    Outputs:
        resampled (channel): (id, hemi: Union['L','R'], resampled: Path) 
            `shape` resampled to `target_sphere` mesh
    */

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

    /*
    Combine two shape files into a surface-based dscalar file

    Arguments:
        id (str): ID key
        left (Path): Left shape file
        right (Path): Right shape file

    Output:
        dscalar (channel): (id, dscalar: Path)
    */


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


workflow resample2native_wf {

    /*
    Perform fsaverage_LR32k to T1w native space resampling using
    MSMSulc-computed registration sphere

    Arguments:
        dscalar (channel): (subject, dscalar: Path) Subject dscalar file in fsaverage_LR32k
        msm_sphere (channel): (subject, hemisphere: Union('L','R'), sphere: Path)
            msm_sphere for a given hemisphere and subject

    Outputs:
        resampled (channel): (subject, resampled: Path) `dscalar` in T1w space
    */
        

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
            split_dscalar.out.left.mix(split_dscalar.out.right)
                              .combine(udscalar.map2sub, by: 0)
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

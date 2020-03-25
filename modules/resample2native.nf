nextflow.preview.dsl = 2

process split_dscalar{

    label 'connectome'

    input:
    tuple val(sub), path(dscalar)

    output:
    tuple val(sub), val('L'), path('L.shape.gii'), emit: left
    tuple val(sub), val('R'), path('R.shape.gii'), emit: right

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
    tuple val(sub), val(hemi), path(shape), path(source_sphere), path(target_sphere)

    output:
    tuple val(sub), val(hemi), path("resampled_shape.${hemi}.shape.gii"), emit: resampled

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
    tuple val(sub), path(left), path(right)

    output:
    tuple val(sub), path("${sub}.dscalar.nii"), emit: dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar \
                !{sub}.dscalar.nii \
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

        //Split dscalar for resampling
        split_dscalar(dscalar)

        //Stack left and right outputs
        left_resample_input = split_dscalar.out.left
                                            .join(msm_sphere, by: [0,1] )
                                            .map{ s,h,d,msm ->    [
                                                                    s,h,d,msm,
                                                                    "${params.atlas}/L.sphere.32k_fs_LR.surf.gii"
                                                                ]
                                                }
        right_resample_input = split_dscalar.out.right
                                            .join(msm_sphere, by: [0,1] )
                                            .map{ s,h,d,msm ->    [
                                                                    s,h,d,msm,
                                                                    "${params.atlas}/R.sphere.32k_fs_LR.surf.gii"
                                                                ]
                                                }
        resample_input = left_resample_input.mix(right_resample_input)
        resample_surf(resample_input)

        //Recombine
        recombine_input = resample_surf.out.resampled
                                        .map { s,h,f -> [s,f] }
                                        .groupTuple ( by: 0, sort: { it.baseName }, size: 2 )
                                        .map { s,f -> [s,f[0],f[1]] }

        recombine(recombine_input)

        emit:
            resampled = recombine.out.dscalar




}

nextflow.preview.dsl = 2

process convert_pial_to_metric{

    label 'connectome'
    input:
    tuple val(sub), val(hemi), path(pial)

    output:
    tuple val(sub), val(hemi),\
    path("${sub}.${hemi}.pial_coordinates.shape.gii"),\
    emit: shape_coords

    shell:
    '''

    wb_command -surface-coordinates-to-metric \
        !{pial} \
        !{sub}.!{hemi}.pial_coordinates.shape.gii

    '''
}

process join_surface_coordinates{

    label 'connectome'
    input:
    tuple val(sub), path(L), path(R)

    output:
    tuple val(sub), path("${sub}.cifti_coords.dtseries.nii"), emit: cifti_coords

    shell:
    '''

    wb_command -cifti-create-dense-scalar \
        -left-metric !{L} -right-metric !{R} \
        !{sub}.cifti_coords.dtseries.nii

    '''
}

// process calculate_centroid{
//
//     label 'connectome'
//     input:
//     tuple val(sub), path(cifti_coords)
//
//     output:
//     tuple val(sub), path("${sub}.centroid.txt"), emit: centroid
//
//     shell:
//     '''
//
//     wb_command -cifti-math \
//         "x * y" \
//         -var "x" !{cifti_coords}
//     '''
//
// }

workflow cortex2scalp_wf{

    take:
        mesh
        roi
        pial

        /*
        The issue is that with the coil optimization placement, a
        ray tracing method will have to be used to identify the nearest cortical region
        for calculating distances

        Meaning that for cortex2scalp distances
        roi_centroids can be used to find the closest scalp distance

        For cortex2scalp distances, you need to use a
        raw tracing approach in order to approximate distances
        */


        /* Pial contains hemisphere information */


        convert_pial_to_metric(pial)
        join_surface_coordinates(
                            .map{s,h,g -> [s,g]}
                            .groupTuple(by: 0, sort: { it.getBaseName() })
                            .map{s,g -> [s, g[0], g[1]]} | view




}

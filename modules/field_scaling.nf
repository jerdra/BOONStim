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
    tuple val(sub), path("${sub}.cifti_coords.dscalar.nii"), emit: cifti_coords

    shell:
    '''

    wb_command -cifti-create-dense-scalar \
        -left-metric !{L} -right-metric !{R} \
        !{sub}.cifti_coords.dscalar.nii

    '''
}

process average_coordinate{

    label 'connectome'

    input:
    tuple val(sub), path(coords), path(roi),\
    val(index), val(ori)

    output:
    tuple val(sub), val(ori), path("${ori}.txt"), emit: avg_coord

    shell:
    '''
    wb_command -cifti-math \
        "x*y" \
        -var "x" !{coords} \
        -select 1 !{index} \
        -var "y" !{roi} \
        !{ori}.dscalar.nii

    wb_command -cifti-roi-average \
        !{ori}.dscalar.nii \
        -cifti-roi !{roi} \
        !{ori}.txt
    '''
}

process make_centroid{

    input:
    tuple val(sub), path(x), path(y), path(z)

    output:
    tuple val(sub), path("${sub}.roi_centroid.txt"), emit: centroid

    shell:
    '''
    paste -d '\n' !{x} !{y} !{z} > !{sub}.roi_centroid.txt
    '''

}

// Calculate ROI --> scalp distance
process calculate_roi2cortex{

    label 'rtms'
    input:
    tuple val(sub), path(mesh), path(centroid)

    output:
    tuple val(sub), path("${sub}.roi_distance.npy"), emit: distance

    shell:
    '''
    /scripts/get_cortex_to_scalp.py !{mesh} --roi !{centroid} !{sub}.roi_distance.npy
    '''
}

process calculate_coil2cortex{

    label 'rtms'
    input:
    tuple val(sub), path(mesh), path(coil_centre)

    output:
    tuple val(sub), path("${sub}.coil_distance.npy"), emit: distance

    shell:
    '''
    /scripts/get_cortex_to_scalp.py !{mesh} --coilcentre !{coil_centre} \
                                    !{sub}.coil_distance.npy
    '''
}

process get_ratio{

    label 'rtms'
    input:
    tuple val(sub), path(cortex2scalp), path(coil2cortex)

    output:
    tuple val(sub), path("${sub}.scaling_factor.txt"), emit: scaling_factor

    shell:
    '''
    #!/usr/bin/env python

    import numpy as np

    c2s = np.load("!{cortex2scalp}")
    c2c = np.load("!{coil2cortex}")
    ratio = c2c/c2s * 100

    with open("!{sub}.scaling_factor.txt","w") as f:
        f.write(str(ratio))
    '''
}

workflow cortex2scalp_wf{

    take:
        mesh
        pial
        roi

    main:
        convert_pial_to_metric(pial)

        // Surface coordinates
        i_join_surface_coordinates = convert_pial_to_metric.out.shape_coords
                                        .map{s,h,g -> [s,g]}
                                        .groupTuple(by: 0, sort: { it.getBaseName() })
                                        .map{s,g -> [s, g[0], g[1]]}
        join_surface_coordinates(i_join_surface_coordinates)

        // Process X,Y,Z separately with the roi
        i_average_coordinate = join_surface_coordinates.out.cifti_coords
                            .join(roi, by: 0)
                            .combine([[1,'X'], [2,'Y'], [3,'Z']])
        average_coordinate(i_average_coordinate)

        // Average the information in each piece
        i_make_centroid = average_coordinate.out.avg_coord
                            .map{ s,o,c -> [s,c] }
                            .groupTuple(by: 0, sort: { it.getBaseName() })
                            .map{ s,f -> [s,f].flatten() }
        make_centroid(i_make_centroid)

        i_calculate_roi2cortex = mesh.join(make_centroid.out.centroid)
        calculate_roi2cortex(i_calculate_roi2cortex)

    emit:
        scalp2cortex = calculate_roi2cortex.out.distance

}

workflow coil2cortex_wf{

    take:
        mesh
        coil_centre

    main:
        i_calculate_coil2cortex = mesh.join(coil_centre)
        calculate_coil2cortex(i_calculate_coil2cortex)

    emit:
        cortex2coil = calculate_coil2cortex.out.distance
}

workflow fieldscaling_wf{

    take:
        mesh
        pial
        roi
        coil_centre

    main:
        cortex2scalp_wf(mesh, pial, roi)
        coil2cortex_wf(mesh, coil_centre)

        i_get_ratio = cortex2scalp_wf.out.scalp2cortex
                                        .join(coil2cortex_wf.out.cortex2coil)
        get_ratio(i_get_ratio)

    emit:
        scaling_factor = get_ratio.out.scaling_factor

}

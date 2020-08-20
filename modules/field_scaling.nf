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

    # Apply the ROI mask to the CIFTI data
    wb_command -cifti-math \
        "x*y" \
        -var "x" !{coords} \
        -select 1 !{index} \
        -var "y" !{roi} \
        !{ori}.dscalar.nii

    # Calculate the average coordinate within the centroid
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

process threshold_roi{

    label 'connectome'

    input:
    tuple val(sub), path(roi)

    output:
    tuple val(sub), path("${sub}.roi_thresholded.dscalar.nii"), emit: mask

    shell:
    '''
    #!/bin/bash

    wb_command -cifti-math "x > 0.5" -var "x" !{roi} \
                            !{sub}.roi_thresholded.dscalar.nii
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

// Alternative formulation
process get_cortical_distance_masked{

    label 'rtms'
    input:
    tuple val(sub), path(mesh), path(left_surf), path(right_surf), path(roi)

    output:
    tuple val(sub), path("${sub}.roi_distance.txt"), emit: distance

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                  --roi !{roi} \
                                  !{sub}.roi_distance.txt
    '''

}

process get_cortical_distance{

    label 'rtms'
    input:
    tuple val(sub), path(mesh), path(left_surf), path(right_surf), path(coilcentre)

    output:
    tuple val(sub), path("${sub}.coil_distance.txt"), emit: distance

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                  --coilcentre !{coilcentre} \
                                  !{sub}.coil_distance.txt
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
                                    !{sub}.coil_distance.txt
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

    c2s = np.genfromtxt("!{cortex2scalp}")
    c2c = np.genfromtxt("!{coil2cortex}")
    ratio = c2c/c2s * 100

    to_write = f"{ratio:.2f}" + "\\n"
    with open("!{sub}.scaling_factor.txt","w") as f:
        f.write(to_write)

    '''
}

process matsimnibs2centre{

    label 'rtms'
    input:
    tuple val(sub), path(matsimnibs)

    output:
    tuple val(sub), path("${sub}_coilcentre.txt"), emit: coil_centre

    shell:
    '''
    #!/usr/bin/env python
    import numpy as np

    matsimnibs = np.load("!{matsimnibs}")
    coil_centre = matsimnibs[:3,-1]
    np.savetxt("!{sub}_coilcentre.txt", coil_centre)
    '''
}

process qc_cortical_distance{

    label 'rtms'
    input:
    tuple val(sub), path(mesh), path(left_surf), path(right_surf),\
    path(coil), path(mask)

    output:
    tuple val(sub), path("${sub}_distqc.geo"), emit: distqc

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                --coilcentre !{coil} --roi !{mask} \
                                --qc-geo /geo/dist.geo \
                                !{sub}_distqc.geo
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0], it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
                }

workflow cortex2scalp_wf{

    take:
        mesh
        pial_left
        pial_right
        roi

    main:

        threshold_roi(roi)
        i_get_cortical_distance = mesh.join(pial_left)
                                      .join(pial_right)
                                      .join(threshold_roi.out.mask)
        get_cortical_distance_masked(i_get_cortical_distance)

    emit:
        scalp2cortex = get_cortical_distance_masked.out.distance

}

workflow coil2cortex_wf{

    take:
        mesh
        pial_left
        pial_right
        coil_centre

    main:

        i_get_cortical_distance = mesh.join(pial_left)
                                      .join(pial_right)
                                      .join(coil_centre)
        get_cortical_distance(i_get_cortical_distance)

    emit:
        cortex2coil = get_cortical_distance.out.distance
}

workflow qc_cortical_distance_wf{
    take:
        mesh
        pial_left
        pial_right
        roi
        coil_centre

    main:

        // Compute mask
        threshold_roi(roi)

        i_qc_cortical_distance = mesh.join(pial_left)
                                 .join(pial_right)
                                 .join(coil_centre)
                                 .join(threshold_roi.out.mask)
        qc_cortical_distance(i_qc_cortical_distance)

}


workflow fieldscaling_wf{

    take:
        mesh
        pial
        roi
        matsimnibs

    main:

        pial.branch(lr_branch).set{pial_surfs}

        // MT calculation
        cortex2scalp_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            roi
        )

        // Coil calculation
        matsimnibs2centre(matsimnibs)
        coil2cortex_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            matsimnibs2centre.out.coil_centre
        )

        i_get_ratio = cortex2scalp_wf.out.scalp2cortex
                                        .join(coil2cortex_wf.out.cortex2coil)
        get_ratio(i_get_ratio)

        qc_cortical_distance_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            roi, matsimnibs2centre.out.coil_centre
        )


    emit:
        scaling_factor = get_ratio.out.scaling_factor

}

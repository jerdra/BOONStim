nextflow.preview.dsl = 2

process convert_pial_to_metric{

    /*
    Convert a surface file into a coordinates shape file
    Arguments:
        sub (str): Subject IDs
        hemi (Union['L','R']): Hemisphere
        pial (Path): .surf.gii pial file

    Outputs:
        shape_coords (channel): (sub, hemi, pial_coords: Path) 
    */

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
    /*
    Join left and right .shape.gii files into a .dscalar.nii file
    Arguments:
        sub (str): Subject ID
        L (path): Path to left .shape.gii 
        R (path): Path to right .shape.gii 

    Outputs:
        cifti_coords (channel): (sub, cifti_coords: Path)
    */

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


process threshold_roi{

    /*
    Apply thresholding to a dscalar file at 0.5
    Arguments:
        uuid (str): Unique ID
        roi (Path): Path to dscalar file to threshold

    Outputs:
        mask (channel): (uuid, thresholded_dscalar: Path)
    */
        

    label 'connectome'

    input:
    tuple val(uuid), path(roi)

    output:
    tuple val(uuid), path("${uuid}.roi_thresholded.dscalar.nii"), emit: mask

    shell:
    '''
    #!/bin/bash

    wb_command -cifti-math "x > 0.5" -var "x" !{roi} \
                            !{uuid}.roi_thresholded.dscalar.nii
    '''
}

process get_cortical_distance_masked{
    /*
    Compute a cortical distance map from an ROI to rest of the brain

    Arguments:
        uuid (str): Unique ID
        mesh (Path): Path to .msh file
        left_surf (Path): Path to left surface GIFTI
        right_surf (Path): Path to right surface GIFTI
        roi (Path): Path to dscalar file

    Outputs:
        distance (channel): (uuid, roi_distance: Path) Distance map dscalar file
    */

    label 'fieldopt'
    input:
    tuple val(uuid), path(mesh), path(left_surf), path(right_surf), path(roi)

    output:
    tuple val(uuid), path("${uuid}.roi_distance.txt"), emit: distance

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                  --roi !{roi} \
                                  !{uuid}.roi_distance.txt
    '''

}

process get_cortical_distance{

    /*
    Compute scalp to cortex distance using a set of ROI coordinates

    Arguments:
        uuid (str): Unique ID
        mesh (Path): Path to .msh file
        left_surf (Path): Path to left GIFTI file
        right_surf (Path): Path to right GIFTI file
        coilcentre (Path): Path to proposed coil location

    Output:
        distance (channel): (uuid, coil_distance: Path) 
    */

    label 'fieldopt'
    input:
    tuple val(uuid), path(mesh), path(left_surf), path(right_surf), path(coilcentre)

    output:
    tuple val(uuid), path("${uuid}.coil_distance.txt"), emit: distance

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                  --coilcentre !{coilcentre} \
                                  !{uuid}.coil_distance.txt
    '''


}

process get_stokes_cf{

    /*
    Compute stokes correction factor from treatment distance and MT distance

    Arguments:
        uuid (str): Unique ID
        cortex2scalp (path): Text file containing MT cortex to scalp distance
        coil2cortex (path): Text file containing treatment coil location to cortex distance

    Outputs:
        stokes_correction (channel): (uuid, stokes_correction: Path) Single value containing stokes correction value
    */

    label 'fieldopt'
    input:
    tuple val(uuid), path(cortex2scalp), path(coil2cortex)

    output:
    tuple val(uuid), path("${uuid}.stokes_correction.txt"), emit: stokes_correction

    shell:
    '''
    #!/usr/bin/env python

    import numpy as np

    c2s = np.genfromtxt("!{cortex2scalp}")
    c2c = np.genfromtxt("!{coil2cortex}")
    cf = 2.8*(c2c - c2s)

    to_write = f"{cf:.2f}" + "\\n"
    with open("!{uuid}.stokes_correction.txt","w") as f:
        f.write(to_write)
    '''
}

process matsimnibs2centre{

    /*
    Extract position coordinates from SimNIBS coil orientation matrix

    Arguments:
        uuid (str): Unique ID
        matsimnibs (Path): Coil orientation matrix

    Outputs:
        coil_centre (channel): (uuid, coil_centre: Path)
    */

    label 'fieldopt'
    input:
    tuple val(uuid), path(matsimnibs)

    output:
    tuple val(uuid), path("${uuid}_coilcentre.txt"), emit: coil_centre

    shell:
    '''
    #!/usr/bin/env python
    import numpy as np

    matsimnibs = np.load("!{matsimnibs}")
    coil_centre = matsimnibs[:3,-1]
    np.savetxt("!{uuid}_coilcentre.txt", coil_centre)
    '''
}

process qc_cortical_distance{

    /*
    Generate QC visualization of distance measurements
    Arguments:
        uuid (str): Unique ID
        mesh (Path): Path to head model .msh file
        left_surf (Path): Path to left surface GIFTI file
        right_surf (Path): Path to right surface GIFTI file
        coil (Path): Coil centre location text file
        mask (Path): dscalar file to display 

    Outputs:
        qc_geo (channel): (uuid, distqc_geo: Path) GMSH file for viewing cortical distance QC
        qc_html (channel): (uuid, distqc_html: Path) Interactive HTML file ffor viewing cortical distance QC
    */

    label 'fieldopt'
    input:
    tuple val(uuid), path(mesh), path(left_surf), path(right_surf),\
    path(coil), path(mask)

    output:
    tuple val(uuid), path("${uuid}_distqc.geo"), emit: qc_geo
    tuple val(uuid), path("${uuid}_distqc.html"), emit: qc_html

    shell:
    '''
    /scripts/cortical_distance.py !{mesh} !{left_surf} !{right_surf} \
                                !{uuid}_distqc.geo \
                                --coilcentre !{coil} --roi !{mask} \
                                --gmsh-qc /geo/dist.geo \
                                --html-qc !{uuid}_distqc.html
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0], it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
                }

workflow cortex2scalp_wf{

    /*
    Compute distance from a cortical region (dscalar ROI) to the scalp using radial measurements
    Arguments:
        mesh (channel): (sub, msh: Path) GMSH files
        pial_left (channel): (sub, left_gii: Path) Left surface GIFTI
        pial_right (channel): (sub, right_gii: Path) Right surface GIFTI
        roi (channel): (sub, dscalar: Path) ROI dscalar file

    Outputs:
        scalp2cortex (channel): (sub, scalp2cortex: Path) Scalp-to-cortex distance
    */
        

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

    /*
    Compute distance from a cortical region (dscalar ROI) to the scalp using radial measurements
    Arguments:
        mesh (channel): (sub, msh: Path) GMSH files
        pial_left (channel): (sub, left_gii: Path) Left surface GIFTI
        pial_right (channel): (sub, right_gii: Path) Right surface GIFTI
        coil_centre (channel): (sub, coil_centre: Path) Coil location

    Outputs:
        coil2cortex (channel): (sub, coil2cortex: Path) Scalp-to-cortex distance
    */

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
    /*
    Generate QC images for scalp to cortex distance measurement

    Arguments:
        mesh (channel): (sub, msh: Path) GMSH files
        pial_left (channel): (sub, left_gii: Path) Left surface GIFTI
        pial_right (channel): (sub, right_gii: Path) Right surface GIFTI
        roi (channel): (sub, dscalar: Path) ROI dscalar file
        coil_centre (channel): (sub, coil_centre: Path) Coil location

    Outputs:
        geo (channel): (uuid, distqc_geo: Path) GMSH file for viewing cortical distance QC
        html (channel): (uuid, distqc_html: Path) Interactive HTML file ffor viewing cortical distance QC
    */
    take:
        mesh
        pial_left
        pial_right
        roi
        coil_centre

    main:

        threshold_roi(roi)

        i_qc_cortical_distance = mesh.join(pial_left)
                                 .join(pial_right)
                                 .join(coil_centre)
                                 .join(threshold_roi.out.mask)
        qc_cortical_distance(i_qc_cortical_distance)

    emit:
        html = qc_cortical_distance.out.qc_html
        geo = qc_cortical_distance.out.qc_geo

}


workflow fieldscaling_wf{
    /*
    Compute stokes correction factor from MT measurements

    Arguments:
        mesh (channel): (sub, msh: Path) GMSH files
        pial (channel): (sub, hemi: Union['L','R'], gii: Path) Pial surface GIFTI file
        roi (channel): (sub, dscalar: Path) ROI dscalar file
        matsimnibs (channel): (sub, matsimnibs: Path) Coil orientation matrix

    Outputs:
        scaling_factor (channel): (uuid, scaling_factor: Path) Stokes correction factor
        geo (channel): (uuid, distqc_geo: Path) GMSH file for viewing cortical distance QC
        html (channel): (uuid, distqc_html: Path) Interactive HTML file ffor viewing cortical distance QC
    */

    take:
        mesh
        pial
        roi
        matsimnibs

    main:

        // Generate unique hash for combination of IDs
        uroi = roi.map{a -> [ a[0], a[0..<-1].join('-').md5(),
                              a[1..<-1], a[-1] ]}
                  .multiMap{s, u, ids, d ->
                    map2id: [u, s, ids]
                    map2dscalar: [u, d]
                    map2sub: [u,s]}

        // Transform inputs to using UUIDs
        pial = uroi.map2sub.map{u,s -> [s,u]}
                .combine(pial, by: 0)
                .map{s,u,h,p -> [u,h,p]}

        mesh = uroi.map2sub.map{u,s -> [s,u]}
                .combine(mesh, by:0)
                .map{s,u,m -> [u,m]}

        matsimnibs = uroi.map2sub.map{u,s -> [s,u]}
                        .combine(matsimnibs, by:0)
                        .map{s,u,m -> [u,m]}

        pial.branch(lr_branch).set{pial_surfs}


        // MT calculation
        cortex2scalp_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            uroi.map2dscalar
        )

        // Coil calculation
        matsimnibs2centre(matsimnibs)
        coil2cortex_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            matsimnibs2centre.out.coil_centre
        )

        i_get_stokes_cf = cortex2scalp_wf.out.scalp2cortex
                                        .join(coil2cortex_wf.out.cortex2coil)
        get_stokes_cf(i_get_stokes_cf)

        // Generate measurement QC pages
        qc_cortical_distance_wf(
            mesh,
            pial_surfs.left,
            pial_surfs.right,
            uroi.map2dscalar, matsimnibs2centre.out.coil_centre
        )

        // Map back to original input
        out_html = uroi.map2id.join(qc_cortical_distance_wf.out.html)
                        .map{u, s, ids, h -> [s, *ids, h]}
        out_geo = uroi.map2id.join(qc_cortical_distance_wf.out.geo)
                        .map{u, s, ids, g -> [s, *ids, g]}

        // Map back into original input
        scaling_output = uroi.map2id.join(get_stokes_cf.out.stokes_correction)
                   .map{u,s,ids,c -> [s,ids,c].flatten()}


    emit:
        scaling_factor = scaling_output
        qc_html = out_html
        qc_geo = out_geo

}

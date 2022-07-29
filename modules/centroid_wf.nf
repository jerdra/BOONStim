// TODO: DEPRECATED REMOVE WORKFLOW??

nextflow.preview.dsl = 2

process split_dscalar {

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

process centroid_project2vol {

    /*
    Project GIFTI shape file into a ribbon-constrained volume

    Arguments:
        sub (str): Subject key
        shape (Path): Input GIFTI file
        pial (Path): Input pial surface
        white (Path): Input white surface
        midthick (Path): Input midthickness surface
        t1 (Path): Target T1 NIFTI

    Output:
        ribbon (channel): (sub, ribbon: Path): Surface data in volume space
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(shape), path(pial), path(white), path(midthick), path(t1)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.ribbon.nii.gz"), emit: ribbon

    shell:
    '''
    wb_command -metric-to-volume-mapping \
                !{shape} \
                !{midthick} \
                !{t1} \
                -ribbon-constrained \
                    !{white} \
                    !{pial} \
                !{sub}.!{hemi}.ribbon.nii.gz
    '''
}

process add_centroid_niftis {

    /*
    Add two NIFTI files together

    Arguments:
        sub (str): Subject ID
        nifti1 (Path): First NIFTI to add
        nifti2 (Path): Second NIFTI to add

    Outputs:
        sumvol (channel): (sub, combined: Path) Summation image
    */

    label 'connectome'

    input:
    tuple val(sub), path(nifti1), path(nifti2)

    output:
    tuple val(sub), path('combined.nii.gz'), emit: sumvol

    shell:
    '''
    wb_command -volume-math \
                "x + y" \
                -var x !{nifti1} \
                -var y !{nifti2} \
                combined.nii.gz
    '''
}

process normalize_vol {

    /*
    Normalize values of the input volume such that they add up to 1

    Arguments:
        sub (str): Subject ID
        vol (Path): NIFTI file to normalize sum on

    Outputs:
        normvol (channel): (subject, normvol: Path) Sum normalized volume NIFTI
    */

    label 'connectome'
    input:
    tuple val(sub), path(vol)

    output:
    tuple val(sub), path('normalized.nii.gz'), emit: normvol

    shell:
    '''
    wb_command -volume-math 'a' fixed.nii.gz -var a !{vol} -fixnan 0

    sum=$(wb_command -volume-stats \
                fixed.nii.gz \
                -reduce SUM)

    wb_command -volume-math \
                "x/${sum}" \
                -var x fixed.nii.gz \
                normalized.nii.gz
    '''

}

process compute_weighted_centroid{

    /*
    Compute centre-of-gravity on an input volume

    Arguments:
        sub (str): Subject ID
        vol (Path): Normalized input volume

    Outputs:
        coord (channel): (subject, coord: Path) Centre of gravity .txt coordinate file
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(vol)

    output:
    tuple val(sub), path("${sub}_ras_coord.txt"), emit: coord

    shell:
    '''
    #!/usr/bin/env python

    import nibabel as nib
    import numpy as np

    #Load image
    img = nib.load("!{vol}")
    affine = img.affine
    data = img.get_data()

    #Mask
    x,y,z = np.where(data > 0)
    coords = np.array([x,y,z])
    vals = data[(x,y,z)]

    #Compute
    weighted_vox = np.dot(coords,vals)[:,np.newaxis]
    r_weighted_vox = np.dot(affine[:3,:3],weighted_vox)
    weighted_coord = r_weighted_vox + affine[:3,3:4]

    #Save
    np.savetxt("!{sub}_ras_coord.txt",weighted_coord)
    '''


}

process get_scalp_seed {

    /*
    Compute a good candidate seed position given a scalar-valued map

    Arguments:
        sub (str): Subject ID 
        mesh (Path): Realistic head model
        dscalar (Path): Scalar-valued map
        l_pial (Path): Left pial GIFTI surface
        r_pial (Path): Right pial GIFTI surface

    Outputs:
        seed (channel): (sub, seed: Path) Candidate coordinate on subject scalp
        qchtml (channel): (sub, html: Path) QC page
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(mesh), path(dscalar), path(l_pial), path(r_pial)

    output:
    tuple val(sub), path("${sub}_seed.txt"), emit: seed
    tuple val(sub), path("${sub}_qcseed.html"), emit: qchtml

    shell:
    '''
    #!/bin/bash

    /scripts/get_scalp_seed.py !{mesh} !{dscalar} !{l_pial} !{r_pial} \
                               !{sub}_seed.txt --qc-file !{sub}_qcseed.html
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0], it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
                }

workflow centroid_radial_wf{

    /*
    Compute approximate radial projection out from pial surface onto scalp
    at a given coordinate

    Arguments:
        msh (channel): (subject, msh: Path) Subject .msh file
        dscalar (channel): (subject, dscalar: Path) dscalar file to compute projection from
        pial (channel): (subject, hemi: Union['L' 'R'], pial: Path) Pial file

    Outputs:
        centroid (channel): (subject, centroid: Path) Computed centroid
        qc (channel): (subject, qchtml: Path) QC file associated with `centroid`
    */
        

    take:
        msh
        dscalar
        pial

    main:

        // Branch out left/right surfaces
        pial.branch(lr_branch).set{pial_surfs}

        i_get_scalp_seed = msh.join(dscalar)
                              .join(pial_surfs.left)
                              .join(pial_surfs.right)

        get_scalp_seed(i_get_scalp_seed)

    emit:
        centroid = get_scalp_seed.out.seed
        qc = get_scalp_seed.out.qchtml


}

workflow centroid_wf{

    /*
    Compute approximate weighted centroid in volume-space from a dscalar file

    Arguments:
        dscalar (channel): (subject, dscalar: Path) dscalar file to compute projection from
        pial (channel): (subject, hemi: Union['L' 'R'], pial: Path) Pial file
        white (channel): (subject, hemi: Union['L' 'R'], white: Path) White file
        midthick (channel): (subject, hemi: Union['L' 'R'], midthickness: Path) Midthickness file
        t1 (channel): (subject, t1: Path) NIFTI file to be used for volume-space projection

    Outputs:
        centroid (channel): (subject, centroid: Path) Computed centroid
    */

    take:
        dscalar
        pial
        white
        midthick
        t1

    main:

        //Split into shapes
        split_dscalar(dscalar)

        //Formulate inputs and mix
        left_project_input = split_dscalar.out.left
                                        .join(pial, by:[0,1])
                                        .join(white, by:[0,1])
                                        .join(midthick, by:[0,1])
                                        .join(t1, by:0)

        right_project_input = split_dscalar.out.right
                                        .join(pial, by:[0,1])
                                        .join(white, by:[0,1])
                                        .join(midthick, by:[0,1])
                                        .join(t1, by:0)

        //Combine into one stream
        project_input = left_project_input.mix(right_project_input)
        centroid_project2vol(project_input)

        //Gather together T1 outputs and sum to form full image
        add_niftis_input = centroid_project2vol.out.ribbon
                                    .groupTuple(by: 0, size: 2)
                                    .map{ s,h,n -> [ s,n[0],n[1] ] }
        add_centroid_niftis(add_niftis_input)

        //Re-normalize
        normalize_vol(add_centroid_niftis.out.sumvol)

        //Calculate centroid
        compute_weighted_centroid(normalize_vol.out.normvol)

    emit:
        centroid = compute_weighted_centroid.out.coord
}

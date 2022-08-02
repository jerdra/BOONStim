nextflow.preview.dsl = 2

process split_dscalar {
    /*
    Split CIFTI dscalar file into constituent
    GIFTI shape files

    Argument:
        sub (str): Subject ID
        dscalar (Path): .dscalar.nii to split

    Outputs:
        left (channel): (sub, 'L', left_shape: Path) Left shape file
        right (channel): (sub, 'R', right_shape: Path) Right shape file
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

process tet_project2vol {
    /*
    Perform ribbon-constrained volume projection

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        shape (Path): Shape file to project to volume
        pial (Path): Pial surface file
        white (Path): White surface file
        midthick (Path): Midthickness surface file
        t1 (Path): Volume T1 to project to

    Outputs:
        ribbon (channel): (sub, hemi, ribbon: Path)
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

process add_tet_niftis {

    /*
    Add NIFTI files together

    Arguments:
        sub (str): Subject ID
        nifti1 (Path): First NIFTI
        nifti2 (Path): Second NIFTI

    Outputs:
        sumvol (channel): (sub, sumvol: Path)
    */
        

    label 'connectome'

    input:
    tuple val(sub), path(nifti1), path(nifti2)

    output:
    tuple val(sub), path("${sub}_combined.nii.gz"), emit: sumvol

    shell:
    '''
    wb_command -volume-math \
                "x + y" \
                -var x !{nifti1} \
                -var y !{nifti2} \
                !{sub}_combined.nii.gz
    '''
}

process tetrahedral_projection {
    /*
    Perform projection from volume space to realistic head model tetrahedrons

    Arguments:
        sub (str): Subject ID
        vol (Path): Volume image to project to realistic head model
        msh (Path): Realistic head model .msh file

    Outputs:
        fem_weights (channel): (sub, fem_weights: Path)
    */
        

    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(vol), path(msh)

    output:
    tuple val(sub), path("${sub}_femfunc.npy"), emit: fem_weights


    shell:
    '''
    /scripts/volume_to_tetrahedral_mapping.py !{vol} !{msh} "!{sub}_femfunc.npy"
    '''


}

workflow tet_project_wf{
    /*

    Tetrahedral projection workflow to map surface-based LR32k data onto
    a realistic head model .msh file

    Arguments:
        dscalar (channel): (sub, dscalar: Path) Input dscalar file to map to realistic head model        
        pial (channel): (sub, pial: Path) Grey matter surface file
        white (channel): (sub, white: Path) White matter surface file
        midthick (channel): (sub, midthick: Path) Midthickness surface file
        t1 (channel): (sub, t1: Path) T1 volume image to use to volume project
        msh (channel): (sub, msh: Path) Realistic head model .msh file

    Outputs:
        fem_weights (channel): (sub, fem_weights: Path) .dscalar.nii file projected into realistic head model

    */

    take:
        dscalar
        pial
        white
        midthick
        t1
        msh

    main:

        split_dscalar(dscalar)

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

        project_input = left_project_input.mix(right_project_input)
        tet_project2vol(project_input)

        add_niftis_input = tet_project2vol.out.ribbon
                                    .groupTuple(by: 0, size: 2)
                                    .map{ s,h,n -> [ s,n[0],n[1] ] }
        add_tet_niftis(add_niftis_input)

        tet_inputs = add_tet_niftis.out.sumvol.join(msh, by: 0)
        tetrahedral_projection(tet_inputs)

    emit:
        fem_weights = tetrahedral_projection.out.fem_weights
}

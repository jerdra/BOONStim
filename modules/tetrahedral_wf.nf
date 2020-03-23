nextflow.preview.dsl = 2

process split_dscalar {

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

    label 'rtms'
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

    get:
        dscalar
        pial
        white
        midthick
        t1
        msh

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
        tet_project2vol(project_input)

        //Gather together T1 outputs and sum to form full image
        add_niftis_input = tet_project2vol.out.ribbon
                                    .groupTuple(by: 0, size: 2)
                                    .map{ s,h,n -> [ s,n[0],n[1] ] }
        add_tet_niftis(add_niftis_input)

        //Tetrahedral projection
        tet_inputs = add_tet_niftis.out.sumvol.join(msh, by: 0)
        tetrahedral_projection(tet_inputs)

        emit:
            fem_weights = tetrahedral_projection.out.fem_weights



}

nextflow.preview.dsl = 2

include { numpy2txt } from "../modules/utils.nf" params(params)

process target_direction {

    /*
    Compute target direction at `coordinate` using Freesurfer .curv files
    
    Arguments:
        sub (str): Subject ID
        coordinate (Path): RAS coordinate file
        fs_dir (Path): Subject Freesurfer Directory

    Outputs:
        direction (channel): (subject, direction: Path)
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(coordinate), path(fs_dir)

    output:
    tuple val(sub), path("${sub}_target_normal.npy"), emit: direction
    tuple val(sub), path("${sub}_target_normal.png"), emit: qc_img
    tuple val(sub), path("${sub}_target_normal.html"), emit: qc_html

    shell:
    """
    python /scripts/get_target_direction.py \
        ${coordinate} \
        ${fs_dir} \
        ${sub}_target_normal.npy \
        --qc-img ${sub}_target_normal.png \
        --qc-html ${sub}_target_normal.html
    """
}

workflow target_direction_wf {
    /*
    Compute radial target direction given a coordinate

    Arguments:
        coordinate (channel):  (subject, coordinate: Path)
        fs_dir (channel): (subject, fs_dir: Path)

    Output:
        direction (channel): (subject, target_direction: Path)
    */
    take:
        coordinate
        fs_dir

    main:
        numpy2txt(coordinate)
        
        target_direction(
            numpy2txt.out.txt
                .join(fs_dir)
        )

    emit:
        direction = target_direction.out.direction
        qc_html = target_direction.out.qc_html
        qc_img = target_direction.out.qc_img
}

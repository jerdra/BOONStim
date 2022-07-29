nextflow.preview.dsl = 2

process target_direction {

    label 'fieldopt'
    label 'bin'

    input:
    tuple val(sub), path(coordinate), path(fs_dir)

    output:
    tuple val(sub), path("${sub}_target_normal.npy"), emit: direction

    shell:
    """
    /scripts/get_target_direction.py \
        ${coordinate} \
        ${fs_dir} \
        ${sub}_target_normal.npy
    """
}

workflow target_direction_wf {
    take:
        coordinate,
        fs_dir

    main:
        target_direction(coordinate.join(fs_dir))

    emit:
        direction = target_direction.out.direction
}

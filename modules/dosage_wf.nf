nextflow.preview.dsl=2

params.use_magnitude = true

process calculate_e100 {

    label 'fieldopt'
    /*
    Calculate E100 metric (100th largest field magnitude)

    Arguments:
        subject (str): Subject ID
        json (Path): Subject target specification JSON
        sim_msh (Path): Simulation mesh
        m2m_dir (Path): M2M directory

    Parameters:
        use_magnitude (bool): Whether to use normE vs direction in json

    Outputs:
        e100 (channel): (subject, e100: Path): E100
    */

    input:
    tuple val(subject), path(json),\
    path(sim_msh), path(m2m_dir)

    output:
    tuple val(subject), path("${subject}_e100.txt"), emit: e100

    script:
    def direction_arg = (params.use_magnitude) ? "" : "--direction-json ${json}"

    """
    /scripts/adjust_dosage.py \
        ${sim_msh} \
        ${m2m_dir} ${direction_arg} \
        ${subject}_e100.txt
    """
}

process scale_didt {
    /*
    Scale didt so that E-field metric in target region
    matches reference

    Arguments:
        subject (str): Subject ID
        metric (Path): Simulation E-field metric
        didt (float): Input dI/dt to compute E-field metric
        reference (float): Reference E-field metric to match

    Outputs:
        didt (channel): (subject, didt: Path) dI/dt to use for subject
    */

    input:
    tuple val(subject), val(didt), path(metric), val(reference)

    output:
    tuple val(subject), path("${subject}_didt.txt"), emit: didt

    shell:
    """
    #!/usr/bin/env python

    import numpy as np
    metric = np.loadtxt('${metric}')
    scaling_factor = reference / np.abs(metric)
    scaled_didt = scaling_factor * ${didt}
    np.savetxt('${subject}_didt.txt', scaled_didt)
    """

}

workflow dosage_adjustment_wf {

    /*
    Workflow to adjust dosage to match a reference E-field value

    Arguments:
        spec_json (channel): (subject: Str, spec: Path) Subject spec json
            must contain: ['dir_X', 'dir_Y', 'dir_Z'] if [use_magnitude=False]
        sim_msh (channel): (subject: Str, msh: Path) Simulation result file
        m2m_dir (channel): (subject: Str, m2m: Path) Subject m2m directory
        reference (value: float): Reference value to match dosage to

    Parameters:
        use_magnitude (value: bool): Whether to match magnitude or direction
            [default=True]
        didt (value: float): dI/dt used in simulations
        max_stim_didt (value: float): Maximum dI/dt of stimulator

    Outputs:
        didt (channel): (subject: Str, didt: Path) Resultant dI/dt value to use
    */

    take:
        spec_sheet
        sim_msh
        m2m_dir
        reference
        use_magnitude

    main:
        adjust_dosage( spec_sheet.join(sim_msh).join(m2m_dir) )
        scale_didt (
            adjust_dosage.out.e100
                .spread([reference])
                .spread([params.didt])
        )

    emit:
        didt = scale_didt.out.didt
}

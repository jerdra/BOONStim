nextflow.preview.dsl=2

process calculate_e100 {

    label 'fieldopt'
    /*
    Calculate E100 metric (100th largest field magnitude)

    Arguments:
        subject (str): Subject ID
        sim_msh (Path): Simulation mesh
        m2m_dir (Path): M2M directory

    Parameters:
        optimize_magnitude (bool): Whether to use normE vs direction in json

    Outputs:
        e100 (channel): (subject, e100: Path): E100
    */

    input:
    tuple val(subject), path(sim_msh), path(m2m_dir)

    output:
    tuple val(subject), path("${subject}_e100.txt"), emit: e100

    shell:
    """
    python /scripts/adjust_dosage.py \
        ${sim_msh} \
        ${m2m_dir} \
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

    label 'fieldopt'

    input:
    tuple val(subject), path(metric), val(reference), val(didt)

    output:
    tuple val(subject), path("${subject}_didt.txt"), emit: didt

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    metric = np.loadtxt('${metric}')
    scaling_factor = ${reference} / np.abs(metric)
    scaled_didt = scaling_factor * ${didt}

    with open('${subject}_didt.txt', 'w') as f:
        f.write(str(scaled_didt))
    """

}

workflow dosage_adjustment_wf {

    /*
    Workflow to adjust dosage to match a reference E-field value

    Arguments:
        sim_msh (channel): (subject: Str, msh: Path) Simulation result file
        m2m_dir (channel): (subject: Str, m2m: Path) Subject m2m directory
        reference (value: float): Reference value to match dosage to

    Parameters:
        optimize_magnitude (value: bool): Whether to match magnitude or direction
            [default=True]
        didt (value: float): dI/dt used in simulations
        max_stim_didt (value: float): Maximum dI/dt of stimulator

    Outputs:
        didt (channel): (subject: Str, didt: Path) Resultant dI/dt value to use
    */

    take:
        sim_msh
        m2m_dir
        reference

    main:
        calculate_e100(sim_msh.join(m2m_dir))
        scale_didt (
            calculate_e100.out.e100
                .spread([reference])
                .spread([params.didt])
        )

    emit:
        didt = scale_didt.out.didt
}

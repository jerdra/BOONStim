nextflow.preview.dsl=2

process adjust_dosage {

    label 'fieldopt'
    label 'bin'
    /*
    Adjust dosage of E-field to match a reference value
    Uses E100 metric (100th largest E-field value)
    
    Arguments:
        subject (str): Subject ID
        json (Path): Subject target specification JSON
        sim_msh (Path): Simulation mesh
        m2m_dir (Path): M2M directory
        reference (float): Reference E-field value (V/m)
        use_magnitude (bool): Whether to use magnitude or not

    Outputs:
        didt (channel): (subject, didt: Path): didt dose to use
    */

    input:
    tuple val(subject), path(json),\
    path(sim_msh), path(m2m_dir), val(reference),\
    val(use_magnitude)

    output:
    tuple val(subject), path("${subject}_didt.txt"), emit: didt

    script:
    def direction_arg = (use_magnitude) ? "" : "--direction-json ${json}"

    """
    /scripts/adjust_dosage.py \
        ${sim_msh} \
        ${m2m_dir} \
        ${reference} ${direction_arg} \
        ${subject}_didt.txt
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
        use_magnitude (value: bool): Whether to match magnitude or direction

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
        adjust_dosage(
            spec_sheet
                .join(sim_msh)
                .join(m2m_dir)
                .spread([reference])
                .spread([use_magnitude])
        )

    emit:
        didt = adjust_dosage.out.didt

            
         


    
}

nextflow.preview.dsl=2

process ciftify_invocation{

    input:
    val sub

    output:
    tuple val(sub), path("${sub}.json"), emit: json

    """

    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${params.ciftify_invocation}'
    x = '${sub}'.replace('sub-','')

    with open(invoke_file,'r') as f:
        j_dict = json.load(f)

    j_dict.update({'participant_label' : [x]})

    with open(out_file,'w') as f:
        json.dump(j_dict,f,indent=4)

    """
}


process ciftify{

    stageInMode 'copy'

    input:
    tuple val(sub), path(json), path("fmriprep/$sub"), path("freesurfer/$sub")

    output:
    tuple val(sub), path("ciftify/${sub}"), emit: ciftify

    shell:
    '''
    mkdir work
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v $(pwd)/work:/work \
    -v !{params.license}:/license \
    -v !{params.resources}:/resources \
    !{params.ciftify_descriptor} $(pwd)/!{json} \
    --imagepath !{params.ciftify} -x --stream
    '''
}

// fMRIPrep anat invocation
process fmriprep_invocation{

    input:
    tuple val(sub), path(invocation)

    output:
    tuple val("$sub"), path("${sub}.json"), emit: json

    shell:
    """

    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${invocation}'
    x = '${sub}'.replace('sub-','')

    print(invoke_file)

    with open(invoke_file,'r') as f:
        j_dict = json.load(f)

    j_dict.update({'participant_label' : [x]})

    with open(out_file,'w') as f:
        json.dump(j_dict,f,indent=4)

    """

}

process fmriprep_anat{

    beforeScript "source /etc/profile"
    module 'slurm'

    input:
    tuple val(sub), path(json)

    output:
    tuple val(sub),\
    path("fmriprep/$sub/ses-01/anat/${sub}_ses-01_*_run-1_desc-preproc_T1w.nii.gz"),\
    emit: preproc_t1


    shell:
    '''
    mkdir work
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v $(pwd)/work:/work \
    -v !{params.license}:/license \
    -v !{params.resources}:/resources \
    !{params.fmriprep_descriptor} $(pwd)/!{json} \
    --imagepath !{params.fmriprep} -x --stream

    '''
}

process run_fmriprep{

    beforeScript "source /etc/profile"
    module 'slurm'

    input:
    tuple val(sub), path(json)

    output:
    tuple val(sub), path("freesurfer/${sub}"), emit: freesurfer
    tuple val(sub), path("fmriprep/${sub}"), emit: fmriprep
    tuple val(sub), path("fmriprep/${sub}.html"), emit: html

    shell:
    '''
    mkdir work
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v $(pwd)/work:/work \
    -v !{params.license}:/license \
    -v !{params.resources}:/resources \
    !{params.fmriprep_descriptor} $(pwd)/!{json} \
    --imagepath !{params.fmriprep} -x --stream

    '''
}

process mri2mesh {

    beforeScript 'source /etc/profile'

    input:
    tuple val(sub), path(t1)

    output:
    tuple val(sub), path("fs_${sub}"), emit: freesurfer
    tuple val(sub), path("m2m_${sub}"), emit: mri2mesh
    tuple val(sub), path("${sub}.geo"), emit: geo
    tuple val(sub), path("${sub}_T1fs_conform.nii.gz"), emit: T1

    // Within container run command
    shell:
    '''
    set +u
    export FS_LICENSE=/license/license.txt
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    source $FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --all !{sub} !{t1}
    '''

}

process update_msh{

    beforeScript 'source /etc/profile'
    label 'gmsh4'

    input:
    tuple val(sub), path("${sub}.geo"), path(m2m)

    output:
    tuple val(sub), path("${sub}.msh"), emit: mesh

    shell:
    '''
    set +u

    sed 's/Merge.*m2m/Merge "m2m/g' !{sub}.geo -i
    /gmsh-sdk/bin/gmsh -3 -format msh2 -o !{sub}.msh !{sub}.geo || true
    '''
}

workflow fmriprep_anat_wf{
    take:
        subs

    main:
        fmriprep_invocation(subs.combine(["$params.fmriprep_anat_invocation"]))
        fmriprep_anat(fmriprep_invocation.out.json)

    emit:
        preproc_t1 = fmriprep_anat.out.preproc_t1
}

 workflow fmriprep_wf{

     take:
         subs

     main:
         fmriprep_invocation(subs.combine(["$params.fmriprep_invocation"]))
         run_fmriprep(fmriprep_invocation.out.json)

     emit:
         fmriprep = run_fmriprep.out.fmriprep
         freesurfer = run_fmriprep.out.freesurfer
         html = run_fmriprep.out.html
 }

// Workflow definition
workflow cifti_meshing_wf {

    take:
        subs

    main:
        // Set up fMRIPrep anatomical pre-processing
        fmriprep_anat_wf(subs)
        mri2mesh(fmriprep_anat_wf.out.preproc_t1)

        // Full fmriprep/ciftify pipeline
        fmriprep_wf(subs)
        ciftify_invocation(subs)
        i_ciftify = ciftify_invocation.out.json
                        .join(fmriprep_wf.out.fmriprep)
                        .join(fmriprep_wf.out.freesurfer)

        // Attach fmriprep outputs
        ciftify(i_ciftify)


        // Add m2m
        update_msh_input = mri2mesh.out.geo.join(mri2mesh.out.mri2mesh)
        update_msh(update_msh_input)

    emit:
        cifti = ciftify.out.ciftify
        freesurfer = fmriprep_wf.out.freesurfer
        fmriprep = fmriprep_wf.out.fmriprep
        fmriprep_html = fmriprep_wf.out.html
        mesh_fs = mri2mesh.out.freesurfer
        mesh_m2m = mri2mesh.out.mri2mesh
        msh = update_msh.out.mesh
        t1fs_conform = mri2mesh.out.T1
}

nextflow.preview.dsl=2

// To allow for optional inputs
params.null_file = "NULL"

def try_fetch_t2 (path) {

    /*
    Attempt to locate files in <path>, select first one if 1 or more findings
    Else return a mock file
    */

    def t2s = file(path)

    def include_t2 = (params.include_t2) ? params.include_t2 : false

    if (!include_t2.toBoolean()) {
        log.info("params.include_t2 set to not true, skipping T2 inclusion")
        file(params.null_file)
    }
    else if (t2s.size() > 1) {
        log.info("Found ${t2s.size()} items for ${path}, selecting first!")
        t2s[0]
    }
    else if (t2s.size() == 0) {
        log.error("No matches for ${path}, skipping T2 inclusion!")
        file(params.null_file)
    }
    else {
        log.info("Found 1 T2 for ${path}")
        t2s[0]
    }
}

process ciftify_invocation{

    /*
    Create subject-specific ciftify invocation file

    Arguments:
        sub (str): Subject ID

    Outputs:
        json (channel): (sub, json: Path) Subject ciftify invocation file
    */

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

    /*
    Subject-specific Ciftify preprocessing pipeline

    Arguments:
        sub (str): Subject ID
        fmriprep (Path): Subject fMRIPrep folder
        freesurfer (Path): Subject Freesurfer folder

    Outputs:
        ciftify (channel): (sub, ciftify: Path) Subject Ciftify folder
        zz_templates (channel): Path ZZ templates folder
        qc_recon (channel): (sub, qc_recon: Path) Subjec recon QC folder
        qc_fmri (channel): (sub, qc_fmri: Path) Subjec fMRI QC folder
    */

    stageInMode 'copy'

    input:
    tuple val(sub), path(json), path("fmriprep/$sub"), path("freesurfer/$sub")

    output:
    tuple val(sub), path("ciftify/${sub}"), emit: ciftify
    path("ciftify/zz_templates"), emit: zz_templates
    tuple val(sub), path("ciftify/qc_recon_all/${sub}"), emit: qc_recon
    tuple val(sub), path("ciftify/qc_fmri/${sub}*"), emit: qc_fmri

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


process fmriprep_invocation{
    /*
    Create subject-specific fmriprep invocation file

    Arguments:
        sub (str): Subject ID

    Outputs:
        json (channel): (sub, json: Path) Subject fmriprep invocation file
    */

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

    /*
    Subject-specific fMRIPrep anatomical only pipeline

    Arguments:
        sub (str): Subject ID
        json (Path): Path to fMRIPrep invocation file with anat_only

    Output:
        preproc_t1 (channel): (sub, preproc: Path) Preprocessed T1 file
    */

    module 'slurm'

    input:
    tuple val(sub), path(json)

    output:
    tuple val(sub),\
    path("${sub}*_desc-preproc_T1w.nii.gz"),\
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
t
    # Find anat file and link to current folder
    find fmriprep/!{sub}/ -type f -name "*preproc_T*w.nii.gz" | \
    grep -v MNI152 | xargs -I [] cp [] .
    '''
}

process run_fmriprep{

    /*
    Subject-specific fMRIPrep pipeline

    Arguments:
        sub (str): Subject ID
        json (Path): Path to fMRIPrep invocation file with anat_only

    Output:
        freesurfer (channel): (subject, fs_dir: Path) Subject freesurfer directory
        fmriprep (channel): (subject, fmriprep: Path) Subject fMRIPrep directory
        html (channel): (subject, fmriprep_html: Path) Subject fMRIPrep HTML file
    */

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

    /*
    Subject-specific SimNIB mri2mesh pipeline

    Arguments:
        sub (str): Subject ID
        t1 (Path): Path to input T1 file
        optional_t2 (Path): Optional T2 input file. [Set to file("NULL") if excluding]

    Output:
        freesurfer (channel): (subject, fs_dir: Path) Subject freesurfer directory
        fmriprep (channel): (subject, fmriprep: Path) Subject fMRIPrep directory
        html (channel): (subject, fmriprep_html: Path) Subject fMRIPrep HTML file
    */

    input:
    tuple val(sub), path(t1), path(optional_t2)

    output:
    tuple val(sub), path("fs_${sub}"), emit: freesurfer
    tuple val(sub), path("m2m_${sub}"), emit: mri2mesh
    tuple val(sub), path("${sub}.geo"), emit: geo
    tuple val(sub), path("${sub}_T1fs_conform.nii.gz"), emit: T1

    script:
    def t2 = (optional_t2.getName() != params.null_file) ? optional_t2 : ""
    """
    set +u
    export FS_LICENSE=/license/license.txt
    source \$FREESURFER_HOME/SetUpFreeSurfer.sh
    source \$FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --all ${sub} ${t1} ${t2}
    """

}

process update_msh{

    /*
    Downgrade the GMSH version of a mesh file to v2
    This is so the file is more usable within the Gmsh python API

    Arguments:
        sub (str): Subject ID
        geo (Path): Subject .geo file
        m2m (Path): Subject mri2mesh m2m directory

    Outputs:
        mesh (channel): (sub, mesh: Path) Path to v2 .msh file
    */


    label 'gmsh4'

    input:
    tuple val(sub), path("${sub}.geo"), path(m2m)

    output:
    tuple val(sub), path("${sub}.msh"), emit: mesh

    shell:
    '''
    set +u

    sed 's/Merge.*m2m/Merge "m2m/g' !{sub}.geo -i
    gmsh -3 -format msh2 -o !{sub}.msh !{sub}.geo || true
    '''
}


process publish_mri2mesh{

    publishDir path: "${params.out}/mri2mesh/${sub}", \
               mode: 'copy', \
               overwrite: true

    input:
    tuple val(sub),\
    path(t1fs), path(m2m), path(fs)

    output:
    tuple val(sub), path(t1fs), path(m2m), path(fs)

    shell:
    '''
    #!/bin/bash
    echo "Transferring !{m2m} and !{fs} to boonstim/!{sub} folder..."
    '''
}

process publish_cifti{

    stageInMode 'copy'

    publishDir path: "$params.out",\
               mode: 'copy'

    input:
    tuple val(sub),\
    path("ciftify/${sub}"), path("ciftify/qc_fmri/*"),\
    path("ciftify/qc_recon_all/${sub}"),\
    path("fmriprep/*"), path(html),\
    path("freesurfer/*"),\
    path("ciftify/zz_templates")

    output:
    tuple path("ciftify/${sub}"),
    path("ciftify/qc_fmri/${sub}*", includeInputs: true),\
    path("ciftify/qc_recon_all/${sub}"),\
    path("fmriprep/${sub}"), path("fmriprep/${sub}.html"),\
    path("freesurfer/${sub}"), path("ciftify/zz_templates")

    shell:
    '''
    echo "Copying fMRIPrep_Ciftify outputs"
    mv !{html} fmriprep/

    '''
}


workflow fmriprep_anat_wf{

    /*
    fMRIPrep anatomical only pipeline

    Arguments:
        subs (channel): Subject IDs

    Parameters:
        bids: Path to BIDS dataset
        fmriprep: Path to fMRIPrep descriptor file
        anat_invocation: Path to fMRIprep anatomical invocation file
        fmriprep_invocation: Path to fMRIprep invocation file
        license: Path to Freesurfer license file
        resources: Path to additional resources folder

    Output:
        preproc_t1 (channel): (sub, preproc: Path) Preprocessed T1 file
    */

    take:
        subs

    main:
        fmriprep_invocation(subs.combine(["$params.anat_invocation"]))
        fmriprep_anat(fmriprep_invocation.out.json)

    emit:
        preproc_t1 = fmriprep_anat.out.preproc_t1
}

 workflow fmriprep_wf{

    /*
    fMRIPrep pipeline

    Arguments:
        subs (channel): Subject IDs

    Parameters:
        bids: Path to BIDS dataset
        fmriprep: Path to fMRIPrep descriptor file
        fmriprep_anat_invocation: Path to fMRIprep anatomical invocation file
        fmriprep_invocation: Path to fMRIprep invocation file
        license: Path to Freesurfer license file
        resources: Path to additional resources folder

    Output:
        freesurfer (channel): (subject, fs_dir: Path) Subject freesurfer directory
        fmriprep (channel): (subject, fmriprep: Path) Subject fMRIPrep directory
        html (channel): (subject, fmriprep_html: Path) Subject fMRIPrep HTML file
    */

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
    /*
    Perform full fMRI data preprocessing and SimNIBS mesh reconstruction (mri2mesh)

    Arguments:
        sub (channel): Subject IDs

    Parameters:
        bids: Path to BIDS dataset
        fmriprep: Path to fMRIPrep descriptor file
        fmriprep_anat_invocation: Path to fMRIprep anatomical invocation file
        fmriprep_invocation: Path to fMRIprep invocation file
        ciftify_descriptor: Path to Ciftify descriptor
        ciftify_invocation: Path to Ciftify invocation file
        license: Path to Freesurfer license file
        resources: Path to additional resources folder
        include_t2 (bool): Include T2w image from BIDS directory in mri2mesh reconstruction

    Outputs:
        cifti (channel): (subject, cifti: Path) Ciftify subject output path
        cifti_qc_fmri (channel): (subject, cifti_qc: Path) Ciftify fmri QC path
        cifti_qc_recon (channel): (subject, cifti_qc: Path) Ciftfy anatomical QC path
        freesurfer (channel): (subject, fs_dir: Path) Subject freesurfer directory
        fmriprep (channel): (subject, fmriprep: Path) Subject fMRIPrep directory
        fmriprep_html (channel): (subject, fmriprep_html: Path) Subject fMRIPrep HTML file
        mesh_fs (channel): (subject, mesh_fs: Path) mri2mesh Freesurfer path
        mesh_m2m (channel): (subject, mesh_m2m: Path) m2m Subject folder
        msh (channel): (subject, mesh_file: Path) Subject .msh file
        t1fs_conform (channel): (subject, t1fs: Path) Subject T1fs conform file
    */

    take:
        subs

    main:
        // Set up fMRIPrep anatomical pre-processing
        fmriprep_anat_wf(subs)

        mri2mesh(
        fmriprep_anat_wf.out.preproc_t1.map { s, t1 ->
            [ s, t1, try_fetch_t2("${params.bids}/${s}/**/*T2w.nii.gz") ]
        })

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

        // Publish outputs
        if (!params.skip_preproc_publish.toBoolean()){
            publish_cifti(
                ciftify.out.ciftify
                    .join(ciftify.out.qc_fmri)
                    .join(ciftify.out.qc_recon)
                    .join(fmriprep_wf.out.fmriprep)
                    .join(fmriprep_wf.out.html)
                    .join(fmriprep_wf.out.freesurfer)
                    .combine(["$params.zz"])
            )

            publish_mri2mesh(
                mri2mesh.out.T1
                    .join(mri2mesh.out.mri2mesh)
                    .join(mri2mesh.out.freesurfer)
            )
        }

    emit:
        cifti = ciftify.out.ciftify
        cifti_qc_fmri = ciftify.out.qc_fmri
        cifti_qc_recon = ciftify.out.qc_recon
        freesurfer = fmriprep_wf.out.freesurfer
        fmriprep = fmriprep_wf.out.fmriprep
        fmriprep_html = fmriprep_wf.out.html
        mesh_fs = mri2mesh.out.freesurfer
        mesh_m2m = mri2mesh.out.mri2mesh
        msh = update_msh.out.mesh
        t1fs_conform = mri2mesh.out.T1
}

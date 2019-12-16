nextflow.preview.dsl = 2

process fmriprep_invocation{

    input:
    val(sub)
        
    output:
    tuple val("$sub"), path("${sub}.json"), emit: json

    shell:
    """
    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${params.anat_invocation}'
    x = '${sub}'.replace('sub-','')

    with open(invoke_file,'r') as f:
        j_dict = json.load(f)

    j_dict.upate({'participant_label' : [x]})

    with open(out_file,'w') as f:
        json.dump(j_dict,f,indent=4)
    """
}

process fmriprep_anat{

    input:
    tuple val(sub), path(json)

    output:
    tuple val(sub), path("fmriprep/${sub}/anat/*preproc_T1w.nii.gz"), emit: preproc_T1

    shell:
    '''
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v !{params.license}:/license \
    !{params.anat_descriptor} $(pwd)/!{json} \
    --imagepath !{params.fmriprep_img} -x --stream
    '''
}

process mri2mesh_brain{

    input:
    tuple val(sub), path(t1)

    output:
    tuple val(sub), path('fs_sub*'), emit: freesurfer
    tuple val(sub), path('m2m_sub*'), emit: m2m
    tuple val(sub), path('sub*_T1fs_conform.nii.gz'), emit: T1

    shell:
    '''
    set +u
    export FS_LICENSE=/license/license.txt
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    source $FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --brain !{sub} !{t1}
    '''
}

process ciftify_invocation{

    input:
    val(sub)
        
    output:
    tuple val("$sub"), path("${sub}.json"), emit: json

    shell:
    """
    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${params.ciftify_invocation}'
    x = '${sub}'.replace('sub-','')

    with open(invoke_file,'r') as f:
        j_dict = json.load(f)

    j_dict.upate({'participant_label' : [x]})

    with open(out_file,'w') as f:
        json.dump(j_dict,f,indent=4)
    """
}

process ciftify{

    input:
    tuple val(sub), path("freesurfer/${sub}"), path(json)

    output:
    tuple val(sub), path("fmriprep/${sub}"), emit: fmriprep
    tuple val(sub), path("ciftify/${sub}"), emit: ciftify
    tuple val(sub), path("freesurfer/${sub}"), emit: freesurfer

    shell:
    '''
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v !{params.license}:/license \
    !{params.ciftify_descriptor} $(pwd)/!{json} \
    --imagepath !{params.ciftify_img} -x --stream
    '''
}

process mri2mesh{

    input:
    tuple val(sub), path(t1), path(t1fs), path(freesurfer), path(m2m)

    output:
    tuple val(sub), path('m2m_sub*'), emit: mri2mesh
    tuple val(sub), path('sub*.geo'), emit: geo

    shell:
    '''
    set +u
    export FS_LICENSE=/license/license.txt
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    source $FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --subcort --head \
            --vol --mni !{sub} !{t1}
    '''
}

process update_msh{

    label 'gmsh4'

    input:
    tuple val(sub), path('sub.geo'), path(m2m)

    output:
    tuple val(sub), path("${sub}.msh"), emit: mesh

    shell:
    '''
    set +u
    sed 's/Merge.*m2m/Merge "m2m/g' sub.geo -i
    /gmsh-sdk/bin/gmsh -3 -bin -format msh2 -o !{sub}.msh sub.geo || true
    '''
}

workflow cifti_meshing {

    get: subs
    main:

        // Step 1: Run fMRIPREP Anat workflow to preprocess and combine data
        fmriprep_invocation(subs)
        fmriprep_anat(fmriprep_invocation.out.json)

        // Pass preproc_T1 into mri2mesh
        mri2mesh_brain(fmriprep_anat.out.preproc_T1)

        // Use freesurfer output from mri2mesh to run ciftify
        ciftify_invocation(subs)
        ciftify_input = mri2mesh_brain.out.freesurfer
                                        .join(ciftify_invocation.out.json, by: 0)
        ciftify(ciftify_input)

        // Continue runs of mri2mesh
        mri2mesh_input = fmriprep_anat.out.preproc_T1
                                        .join(mri2mesh_brain.out.T1, by: 0)
                                        .join(mri2mesh_brain.out.freesurfer, by: 0)
                                        .join(mri2mesh_brain.out.m2m, by: 0)
        mri2mesh(mri2mesh_input)
        update_msh(mri2mesh.out.geo)

    emit:
       cifti = ciftify.out.ciftify 
       fmriprep = ciftify.out.fmriprep
       mesh_fs = mri2mesh_brain.out.freesurfer
       mesh_m2m = mri2mesh.out.m2m
       msh = update_msh.out.mesh
       t1fs_conform = mri2mesh_brain.out.T1
}


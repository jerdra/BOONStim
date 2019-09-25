nextflow.preview.dsl=2


// Process definitions
process fmriprep_invocation{

    input:
    val sub
    
    output:
    tuple val("$sub"), path("${sub}.json"), emit: json

    echo "true"
    
    """

    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${params.anat_invocation}'
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
    tuple val(sub), path("fmriprep/$sub/anat/*preproc_T1w.nii.gz"), emit: preproc_t1
    tuple val(sub), path("fmriprep"), emit: fmriprep

    shell:
    '''
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v !{params.license}:/license \
    !{params.anat_descriptor} $(pwd)/!{json} \
    --imagepath !{params.simg} -x --stream

    #Make available to execution working directory to run MRI2MESH
    cp !{params.out}/fmriprep/!{sub}/anat/*preproc_T1w.nii.gz $(pwd)

    '''
}


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
    
    input:
    tuple val(sub), path(json), path(fmriprep)

    output:
    tuple val(sub), path('fmriprep'), emit: fmriprep
    tuple val(sub), path('ciftify'), emit: ciftify
    tuple val(sub), path('freesurfer'), emit: freesurfer

    shell:
    '''
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v $(pwd):/output \
    -v !{params.license}:/license \
    !{params.ciftify_descriptor} $(pwd)/!{json} \
    --imagepath !{params.simg} -x --stream
    '''

    

}

process mri2mesh {
    
    beforeScript 'source /etc/profile'
    
    input:
    tuple val(sub), path(t1)

    output:
    tuple val(sub), path('fs_sub*'), emit: freesurfer
    tuple val(sub), path('m2m_sub*'), emit: mri2mesh
    tuple val(sub), path('sub*.geo'), emit: geo

    // Within container run command
    shell:
    '''
    tuple +u
    export FS_LICENSE=/license/license.txt
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    source $FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --all !{sub} !{t1}
    '''

}

process update_msh{

    beforeScript 'source /etc/profile'
    
    input:
    tuple val(sub), path(geo)
    
    output:
    path 'sub.msh', emit: mesh

    shell:
    '''
    tuple +u
    
    /gmsh-sdk/bin/gmsh -3 -bin -format msh2 -o sub.msh sub.geo
    '''
}

// Workflow definition
workflow cifti_meshing {

    //Subject list as inputs with implicit input/output dir params
    get: subs 
    main:
    fmriprep_invocation(subs)
    ciftify_invocation(subs)

    fmriprep_anat(fmriprep_invocation.out.json)[1].view()

    //// Join
    //ciftify_inp = fmriprep_anat.out.fmriprep

    //ciftify(ciftify_invocation.out.json)

    //mri2mesh(fmriprep_anat.out.preproc_t1)
    //update_msh(mri2mesh.out.geo)

    ////Organize channels for output

}

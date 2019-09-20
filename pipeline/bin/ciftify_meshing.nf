// BOONStim Nextflow Pipeline implementation

//Main argument is BIDS directory and subject
if (!params.bids || !params.out){

    println('Insufficient input specification!')
    println('Needs --bids, --out!')
    println('Exiting...')
    System.exit(1)

}

// Enable if using T2 images
params.use_t2 = false

println("BIDS Directory: $params.bids")
println("Output Directory: $params.out")

//Now check BOSH invocation of Ciftify
if (!params.anat_invocation || !params.ciftify_invocation || !params.ciftify_descriptor || !params.anat_descriptor) {

    println('Missing BOSH invocations and descriptor for fmriprep-ciftify!')
    println('Please have an --anat-only invocation for fmriprep and an associated descriptor')
    println('Exiting with Error')
    System.exit(1)

}


println("Using Descriptor Files: $params.anat_descriptor and $params.ciftify_descriptor")
println("Using Invocation Files: $params.anat_invocation and $params.ciftify_invocation")

///////////////////////////////////////////////////////////////////////////////

// Main Processes

all_dirs = file(params.bids)
bids_channel = Channel
                    .from(all_dirs.list())
                    .filter { it.contains('sub') }
                    .take(1)


process modify_invocations_for_anat{
    
    // Takes a BIDS subject identifier
    // Modifies the template invocation json and outputs
    // subject specific invocation

    input:
    val sub from bids_channel

    output:
    file "${sub}.json" into invoke_anat_json
    val "${sub}" into invoked_anat_sub

    """

    #!/usr/bin/env python

    import json
    import sys

    out_file = '${sub}.json'
    invoke_file = '${params.anat_invocation}'
    x = '${sub}'.replace('sub-','')

    with open(invoke_file,'r') as f:
        j_dict = json.load(f)
    
    j_dict.update({'participant_label' : [x]})

    with open(out_file,'w') as f:
        json.dump(j_dict,f,indent=4)

    """ 
}

//FMRIPREP anatomical workflow
process run_anat_fmriprep{

    beforeScript "source /etc/profile"
    scratch true
    module 'slurm'
    echo true
    
    input:
    file sub_input from invoke_anat_json
    val sub from invoked_anat_sub

    output:
    file "*preproc_T1w.nii.gz" into preproc_T1
    val "${sub}" into complete_anat_sub
    

    shell:
    '''

    bosh exec launch \
    -v !{params.bids}:/bids \
    -v !{params.out}:/output \
    -v !{params.license}:/license \
    !{params.anat_descriptor} $(pwd)/!{sub_input} \
    --imagepath !{params.simg} -x --stream

    #Make available to execution working directory to run MRI2MESH
    cp !{params.out}/fmriprep/!{sub}/anat/*preproc_T1w.nii.gz $(pwd)

    '''

}

complete_anat_sub.into{ input_ciftify_sub; input_mri2mesh_sub }

// Make invocation for ciftify
process modify_invocations_for_ciftify{


    input:
    val sub from input_ciftify_sub

    output:
    file "${sub}.json" into invoke_ciftify_json
    val "${sub}" into invoked_ciftify_sub

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

// Continue run using Ciftify
process run_ciftify{

    beforeScript "source /etc/profile"
    scratch true
    module 'slurm'
    echo true

    input:
    file sub_input from invoke_ciftify_json

    shell:
    '''
    bosh exec launch \
    -v !{params.bids}:/bids \
    -v !{params.out}:/output \
    -v !{params.license}:/license \
    !{params.ciftify_descriptor} $(pwd)/!{sub_input} \
    --imagepath !{params.simg} -x --stream
    '''

}

// Run MRI2MESH command
process run_mri2mesh{

    beforeScript "source /etc/profile"
    module 'slurm'
    echo true

    publishDir "${params.out}/sim_mesh/${sub}/", mode: 'move', \
                pattern: "sub*!(.geo)"

//    publishDir "${params.out}/sim_mesh/${sub}/", mode: 'copy', \
//                pattern: "sub*.geo"

    containerOptions "-B ${params.license}:/license"

    input: 
    file t1 from preproc_T1
    val sub from input_mri2mesh_sub
    
    //Define output
    output:
    file 'fs_sub*' into fs_sub
    file 'm2m_sub*' into m2m_sub
    set val(sub), file('sub*.geo') into mesh_files
    
    // Within container run command
    shell:
    '''
    set +u
    export FS_LICENSE=/license/license.txt
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    source $FSLDIR/etc/fslconf/fsl.sh
    mri2mesh --all !{sub} !{t1}
    rm !{t1}

    #GMSH3 .msh is bad
    rm !{sub}.msh
    '''
}


// Use newest GMSH to formulate a msh2 file
process update_msh{

    beforeScript "source /etc/profile"
    echo true
    publishDir "${params.out}/sim_mesh/${sub}/${sub}.msh", mode: 'move'

    input:
    set val(sub), file("sub.geo") from mesh_files

    output:
    file 'sub.msh' into updated_mesh

    shell:
    '''
    set +u
    
    /gmsh-sdk/bin/gmsh -3 -bin -format msh2 -o sub.msh sub.geo
    
    '''

}

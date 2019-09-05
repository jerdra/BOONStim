
// Prepares the weight function on the mesh surface and the domain for optimization.


// INPUTS
// out                      Output folder containing ciftify and meshing outputs

// IMPLICIT configuration variables (overrideable)
// clean_cfg                Surface-based cleaning configuration file

// OUTPUTS
// Numpy file with quadratic surface constants
// Numpy file containing nodal weights
// Transformation affine for mapping quadratic domain inputs to sampling surface
// Nodal weights

// INPUT SPECIFICATION
//Check input parameters
if (!params.out){

    log.info('Insufficient input specification!')
    log.info('Needs --out!')
    log.info('Exiting...')
    System.exit(1)

}

//if (!params.template_dir) {
//
//    log.info("Insufficient input specification!"    )
//    log.info('Needs --template_dir!')
//    log.info('Exiting...')
//    System.exit(1)
//
//}
//////////////////////////////////////////////////////////

// MAIN PROCESSES

// Get a list of all subjects with mesh generated
println(params.out) 
sim_mesh_dirs = "$params.out/ciftify/sub-*"

weightfunc_input = Channel.create()
affine_extract_input = Channel.create()

weightfunc_subs = Channel.fromPath(sim_mesh_dirs, type: 'dir')
                    .map { n -> n.getBaseName() }
                    .tap ( weightfunc_input )
                    .map { n -> [
                                    n,
                                    file("$params.out/sim_mesh/$n/${n}_T1fs_conform.nii.gz")
                                ]
                         }
                    .tap ( affine_extract_input )

// Compute weighting function to use give the output file!
process compute_weightfunc {


    module 'connectome-workbench'
    module 'python/3.6.3-anaconda-5.0.1'

    input:
    val sub from weightfunc_input

    output:
    set val(sub), file("*${params.outfile}") into weightfunc_outputs

    //Run custom function
    shell:
    '''
    !{params.weightfunc} !{params.out} !{sub} !{params.outfile}
    '''

}

// If more than one file, take the average between them
process combine_weightfiles {

    echo true
    module 'connectome-workbench'
    module 'python/3.6.3-anaconda-5.0.1'

    input:
    set val(sub), file("input*") from weightfunc_outputs

    output:
    set val(sub), file("combined_${params.outfile}") into weightfiles

    
    shell:
    '''
    find input* | xargs -I {} echo -cifti {} | \
       xargs wb_command -cifti-average combined_!{params.outfile}
    '''

}

//Resample data into SimNIBS space
process weightfunc_to_tetra {

    echo true
    module 'connectome-workbench'

    input:
    set val(sub), file("combined_${params.outfile}") from weightfiles
    file output from file(params.out)

    output:
    set val(sub), file("weight.mesh.dscalar.nii") into resamp_weightfiles

    shell:
    '''
    #Split up dscalar file (which contains volume but like, idk rn)
    wb_command -cifti-separate \
                combined_!{params.outfile} \
                COLUMN \
                -metric CORTEX_RIGHT weight.R.shape.gii \
                -metric CORTEX_LEFT  weight.L.shape.gii

    #Do resampling
    mninonlin=!{output}/ciftify/!{sub}/MNINonLinear/fsaverage_LR32k/
    registration=!{output}/registration/!{sub}/
    wb_command -metric-resample \
                weight.R.shape.gii \
                $mninonlin/!{sub}.R.sphere.32k_fs_LR.surf.gii \
                $registration/!{sub}.R.sphere.reg_msm.surf.gii \
                BARYCENTRIC \
                weight.R.mesh.shape.gii 

    wb_command -metric-resample \
                weight.L.shape.gii \
                $mninonlin/!{sub}.L.sphere.32k_fs_LR.surf.gii \
                $registration/!{sub}.L.sphere.reg_msm.surf.gii \
                BARYCENTRIC \
                weight.L.mesh.shape.gii


    #Combine and pass? 
    wb_command -cifti-create-dense-scalar    \
                -left-metric weight.L.mesh.shape.gii \
                -right-metric weight.R.mesh.shape.gii \
                weight.mesh.dscalar.nii
                
    '''

}

com_inputs = Channel.create()
tetra_inputs = Channel.create()
resamp_weightfiles.into { com_inputs; ribbon_input}


// Using the weight function, find the centre of mass using a user-provided function
process calculate_CoM {

    echo true
    module 'connectome-workbench'
    module 'python/3.6.3-anaconda-5.0.1'

    input:
    set val(sub), file("weight.mesh.dscalar.nii") from com_inputs
    file output from file(params.out)
    
    output:
    set val(sub), file("centre_coordinate.txt") into com_out


    shell:
    '''
    !{params.massfunc} weight.mesh.dscalar.nii !{output} !{sub} "centre_coordinate.txt"
    '''

}


// Extract affine transformation matrix using ciftify output
process extract_affine {

    input:
    set val(sub), file("t1fs.nii.gz") from affine_extract_input

    output:
    set val(sub), file("affine.npy") into affines
    

    shell:
    '''
    #!/usr/bin/env python


    import nibabel as nib
    import numpy as np
    import os

    img = nib.load("t1fs.nii.gz")
    affine = img.affine    
    np.save("affine.npy",affine)
    '''
}


process ribbon_projection { 
    
    module 'connectome-workbench'
    
    input:
    set val(sub), file("weight.dscalar.nii") from ribbon_input
    file "output" from file(params.out)

    output:
    set val(sub), file("ribbon.nii.gz") into ribbon_out

    shell:
    '''
    #!/bin/bash
    
    #Split CIFTI file into GIFTI
    wb_command -cifti-separate \
                weight.dscalar.nii \
                COLUMN \
                -metric CORTEX_RIGHT weight.R.shape.gii \
                -metric CORTEX_LEFT weight.L.shape.gii

    #Project each half into volume space
    wb_command -metric-to-volume-mapping \
                weight.R.shape.gii \
                output/registration/!{sub}/!{sub}.R.midthickness.surf.gii \
                output/sim_mesh/!{sub}/!{sub}_T1fs_conform.nii.gz \
                -ribbon-constrained \
                output/registration/!{sub}/!{sub}.R.white.surf.gii \
                output/registration/!{sub}/!{sub}.R.pial.surf.gii \
                ribbon.R.nii.gz

    wb_command -metric-to-volume-mapping \
                weight.L.shape.gii \
                output/registration/!{sub}/!{sub}.L.midthickness.surf.gii \
                output/sim_mesh/!{sub}/!{sub}_T1fs_conform.nii.gz \
                -ribbon-constrained \
                output/registration/!{sub}/!{sub}.L.white.surf.gii \
                output/registration/!{sub}/!{sub}.L.pial.surf.gii \
                ribbon.L.nii.gz

    #Squish them together
    wb_command -volume-math \
                "a+b" \
                -var a ribbon.R.nii.gz \
                -var b ribbon.L.nii.gz \
                ribbon.nii.gz
    '''


}

tetra_input = ribbon_out
                    .map { n,r ->   [
                                        n,
                                        r,
                                        file("$params.out/sim_mesh/$n/${n}.msh")
                                    ]
                         }

//Should output data.msh, tetraweights
process tetrahedral_projection {

    publishDir "$params.out/fem_optimization/$sub/", \
                saveAs: { "${sub}_$it" }, \
                mode: 'copy'

    input:
    set val(sub), file("ribbon.nii.gz"), file("data.msh") from tetra_input

    output:
    set val(sub), file("tetraweights.npy") into tetra_out

    shell:
    """
    $params.rtms_bin/volume_to_tetrahedral_mapping.py ribbon.nii.gz data.msh tetraweights.npy
    """

}

// Make a channel with mesh, affine, and centroid inputs
surface_patch_input = com_out
                        .map { s,c ->   [
                                            s,
                                            c,
                                            file("$params.out/sim_mesh/$s/${s}.msh")
                                        ]
                             }
                        .join ( affines )
//                        .subscribe { log.info("$it") }


// Make a surface patch 
process make_surface_patch {

    input:
    set val(sub), file("coordinate.txt"), file("data.msh"), file("affine.npy") from surface_patch_input

    output:
    set val(sub), file("patch_dilated_coords.npy"), file("patch_mean_norm.npy") into surface_patch

    """
    $params.rtms_bin/extract_surface_patch.py "data.msh" "affine.npy" "coordinate.txt" "patch" 
    """

}

// Parameterize the surface patch using quadratic fit
//
process parameterize_surface {

    publishDir "$params.out/fem_optimization/$sub/", \
                saveAs: { "${sub}_$it" }, \
                mode: 'copy'

    input:
    set val(sub), file("patch.npy"), file("norm.npy") from surface_patch

    output:
    set val(sub), file("surf_C.npy"), file("surf_R.npy"), file("surf_bounds.npy") into param_surf

    """
    $params.rtms_bin/parameterize_surface_patch.py "patch.npy" "norm.npy" "surf"
    """

}

//FINAL RELEVANT OUTPUTS FOR OPTIMIZATION....
// PARAM_SURF
// TETRA_OUT

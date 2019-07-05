#!/bin/bash

// Registration script from Freesurfer --> Ciftify


//INPUTS:
//out                   Directory containing outputs for sim_mesh
//bids                  Directory with BIDS subjects

//IMPLICIT CONFIG VARIABLES:
//template_dir          Directory containing templates for atlases


//OUTPUTS:
//Script will output into output_dir/sim_mesh/sub/registration


//Check input parameters
if (!params.out){

    println('Insufficient input specification!')
    println('Needs --out!')
    println('Exiting...')
    System.exit(1)

}

///////////////////////////////////////////////////////////////////////////////

// MAIN PROCESSES

// Pull available output directories
sim_mesh_dir = params.out + '/sim_mesh/sub-*'
sub_channel = Channel
                .fromPath(sim_mesh_dir, type:'dir' ) 

//Need a list with subject directories and an L/R version of each
sphere_tuple = sub_channel
                    .map { n -> [ n.baseName, n ] }
                    .map { n -> [n[0], n[1] + "/fs_${n[0]}/surf/"]}
                    .spread( ['L','R'] )
                    .spread( ['sphere', 'sphere.reg'] )
                    .map { n -> [ 
                        n[0],
                        n[1] + "/${n[2].toLowerCase()}h.${n[3]}",
                        n[2] + '.' +n[3]
                        ] }

//Convert from freesurfer --> gifti
process convert_spheres2gifti {

    label 'freesurfer'
    stageInMode 'copy'


    containerOptions "-B ${params.license}:/license"

    input:
    set val(sub), file(sphere), val(output) from sphere_tuple

    output:
    file "${output}.surf.gii" into fs_derived_sphere_files
    val sub into fs_derived_sphere_subject

    shell:
    '''
    set +u 

    echo !{output}
    export FS_LICENSE=/license/license.txt
    mris_convert !{sphere} !{sphere}.surf.gii        
    mv !{sphere}.surf.gii !{output}.surf.gii
    '''

}

//Pull structure from output [ L/R, file object ]
structure_map = ['L': 'CORTEX_LEFT',
                 'R': 'CORTEX_RIGHT']
test = fs_derived_sphere_files
                .map { n -> [ structure_map[n.baseName.take(1)], n]}
                .merge ( fs_derived_sphere_subject )


//Assign connectome workbench names
process assign_surface_properties {

    label 'connectome'
    echo true
    stageInMode 'copy'

    publishDir "${params.out}/registration/$sub/", \
                mode: 'copy', \
                saveAs: { "${sub}.${sphere}" }
    
    input:
    set val(hemi), file(sphere), val(sub)  from test

    output:
    file "$sphere" into assigned_spheres

    """
    wb_command -set-structure ${sphere} ${hemi} -surface-type "SPHERICAL"
    """

}

//Spherical Deformation method


// Registration script from Freesurfer --> Ciftify


//INPUTS:
//out                   Directory containing outputs for sim_mesh

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

//Associate 

//OUTPUT: [sub, path2sphere, output_sphere]
sphere_tuple = sub_channel
                .map { n -> [ n.baseName, n ] }
                .map { n -> [n[0], n[1] + "/fs_${n[0]}/surf/"]}
                .spread( ['L','R'] )
                .tap { sulcal }
                .spread( ['sphere', 'sphere.reg'] )
                .map { n -> [ 
                    n[0],
                    n[1] + "/${n[2].toLowerCase()}h.${n[3]}",
                    n[2] + '.' +n[3]
                    ] }

//Extract sulcal and white
//OUTPUT: [sub, hemi.sulc, hemi.white, hemi]
sulcal = sulcal
            .map { n -> [ n[0],
                          n[2],
                          n[1] + "/${n[2].toLowerCase()}h.sulc",
                          n[1] + "/${n[2].toLowerCase()}h.white"
                        ]
                }

structure_map = ['L': 'CORTEX_LEFT',
                 'R': 'CORTEX_RIGHT']


//Convert sulcal height map into metric
process convert_sulcal {

    label 'freesurfer'
    containerOptions "-B ${params.license}:/license" 

    input:
    set val(sub), val(hemi), file(sulc), file(white)  from sulcal

    output:
    set val(sub), val(hemi), file("${hemi}.sulc.native.shape.gii") into sulcal_metric


    """
    set +u
    export FS_LICENSE=/license/license.txt
    mris_convert -c $sulc $white ${hemi}.sulc.native.shape.gii
    """


}

//Add structure info using hemisphere label
sulcal_struct_input = sulcal_metric
                            .map { n -> [ 
                                         n[0],
                                         n[1],
                                         structure_map[n[1]], n[2] 
                                        ]
                                }


//Transform sulcal height --> depth, and add structure
process invert_sulcal_add_struct {

    label 'connectome'

    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.$it" }, \
                mode: 'copy'

    input:
    set val(sub), val(hemi), val(structure), file(sulc) from sulcal_struct_input

    output:
    set val(sub), val(hemi), file(sulc) into sulcal_depth

    """
    wb_command -metric-math 'a*(-1)' -var 'a' $sulc $sulc
    wb_command -set-structure $sulc $structure
    """

}
process convert_spheres2gifti {

    label 'freesurfer'
    stageInMode 'copy'

                
    containerOptions "-B ${params.license}:/license" 

    input:
    set val(sub), file(sphere), val(output) from sphere_tuple

    output:
    //file "${output}.surf.gii" into fs_derived_sphere_files
    //val sub into fs_derived_sphere_subject
    set val(sub), file("${output}.surf.gii") into fs_derived_spheres

    shell:
    '''
    set +u 

    export FS_LICENSE=/license/license.txt
    mris_convert !{sphere} !{sphere}.surf.gii        
    mv !{sphere}.surf.gii !{output}.surf.gii
    '''

}



//OUTPUT: [sub, surf.gii, structure]
spheres_2_assign = fs_derived_spheres
                        .map{ n -> [ n[0],
                                     n[1],
                                     structure_map[n[1].name.take(1)] 
                                   ] 
                            }

process assign_surface_properties {

    label 'connectome'
    stageInMode 'copy'

    publishDir "${params.out}/sim_mesh/$sub/registration/", \
                mode: 'copy', \
                saveAs: { "${sub}.${sphere}" }
    
    input:
    set val(sub), file(sphere), val(structure) from spheres_2_assign

    output:
    set val(sub), file(sphere) into assigned_spheres

    """
    wb_command -set-structure ${sphere} ${structure} -surface-type "SPHERICAL"
    """

}


reg_spheres = Channel.create()
native_spheres = Channel.create()

//Split into 2 channels based on reg vs non-reg spheres
assigned_spheres
        .choice(reg_spheres,
                native_spheres)
               { a -> a[1].name.contains('reg') ? 0 : 1 }

//OUTPUT: [sub, reg_sphere, hemi]
reg_spheres = reg_spheres
                    .map { n -> [
                                    n[0],
                                    n[1],
                                    n[1].name.take(1)
                                ]
                        }

native_spheres = native_spheres
                    .map { n -> [
                                    n[0],
                                    n[1],
                                    n[1].name.take(1)
                                ]
                        }


// Spherical Deformation method
process spherical_deformation {

    label 'connectome'
    stageInMode 'copy'
    
    containerOptions "-B ${params.atlasdir}:/atlas"

    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.$it" }, \
                mode: 'copy'
    input:
    set val(sub), file(reg_sphere), val(hemi) from reg_spheres

    output:
    set val(sub), file("${hemi}.sphere.reg.reg_LR.native.surf.gii"), val(hemi) into reg_LR_spheres
    

    """
    wb_command -surface-sphere-project-unproject \
    ${reg_sphere} \
    /atlas/fsaverage.${hemi}.sphere.164k_fs_${hemi}.surf.gii \
    /atlas/fs_${hemi}-to-fs_LR_fsaverage.${hemi}_LR.spherical_std.164k_fs_${hemi}.surf.gii \
    ${hemi}.sphere.reg.reg_LR.native.surf.gii
    """ 
}


//Merge sphere.reg.reg_LR and native sphere by subject and hemisphere
//OUTPUT: [sub, hemi, native_sphere, reg_LR_sphere
native_and_reg_spheres = native_spheres
                            .join(reg_LR_spheres, by : [0,2])

//Perform affine regression for surface matching
process spherical_affine_regression {

    label 'connectome'
    stageInMode 'copy'

    input:
    set val(sub), val(hemi), file(sphere), file(reg_LR_sphere) from native_and_reg_spheres

    output:
    set val(sub), val(hemi), file(sphere), file("${hemi}_affine.mat") into affine_spheres

    
    """
    wb_command -surface-affine-regression \
    ${sphere} \
    ${reg_LR_sphere} \
    ${hemi}_affine.mat
    """
    

}

affine_spheres
        .tap { affines }
        .map { n -> [ n[0], n[1], n[2] ] }
        .tap { spheres2rot }

affines = affines
                .map { n -> [ n[0], n[1], n[3] ] }

//Perform normalization of rotation
process normalize_rotation {

    label 'optimize'
    stageInMode 'copy'

    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.${hemi}.affine.mat" }, \
                mode: 'copy'

    input:
    set val(sub), val(hemi), file(affine) from affines

    output:
    set val(sub), val(hemi), file("norm_affine.mat") into normalized_affines

    shell:
    '''
    #!/usr/bin/env python
    import numpy as np

    M = np.genfromtxt("!{affine}")
    M[:,3] = 0
    M[3,3] = 1
    
    linear_map = M[:3,:3]
    U,S,V = np.linalg.svd(linear_map)
    M[:3,:3] = np.matmul(U,V)
    np.savetxt("norm_affine.mat",M)
    
    '''
}

//Join back affines --> native spheres
rotation_inputs = spheres2rot
                        .join(normalized_affines, by: [0,1])

//Apply affine transform to native spheres
process apply_affine {

    label 'connectome'
    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.$it" }, \
                mode: 'copy'

    input:
    set val(sub), val(hemi), file(sphere), file(affine) from rotation_inputs

    output:
    set val(sub), val(hemi), file("${hemi}.sphere_rot.surf.gii") into rotated_sphere
   
    """
    wb_command -surface-apply-affine \
    $sphere $affine ${hemi}.sphere_rot.surf.gii

    wb_command -surface-modify-sphere \
    ${hemi}.sphere_rot.surf.gii 100 ${hemi}.sphere_rot.surf.gii
    """

}



//Perform MSMSulc registration!

//Need:
//rot sphere
//native sulc
//reference sphere
//reference sulc

msm_inputs =  rotated_sphere
                        .join(sulcal_depth, by: [0,1])
                        .map { n -> [
                                        n[0],
                                        n[1],
                                        structure_map[n[1]],
                                        n[2],
                                        n[3]
                                    ]
                            }

process msm_sulc {

    label 'connectome'
    stageInMode 'copy'
    containerOptions "-B ${params.atlasdir}:/atlas -B ${params.msm_conf}:/msm_conf"

    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.$it" }, \
                mode: 'copy'

    input:
    set val(sub), val(hemi), val(structure), file(sphere), file(sulc) from msm_inputs
    
    output:
    set val(sub), val(hemi), file(sphere), file("${hemi}.sphere.reg_msm.surf.gii") into msm_outputs

    """
    /msm/msm --inmesh=$sphere \
        --indata=$sulc \
        --refmesh=/atlas/fsaverage.${hemi}_LR.spherical_std.164k_fs_LR.surf.gii \
        --refdata=/atlas/${hemi}.refsulc.164k_fs_LR.shape.gii \
        --conf=/msm_conf/MSMSulcStrainFinalconf \
        --out=${hemi}. \
        --verbose

    mv "${hemi}.sphere.reg.surf.gii" \
        "${hemi}.sphere.reg_msm.surf.gii"

    wb_command -set-structure ${hemi}.sphere.reg_msm.surf.gii \
                                ${structure}
    
    """
}


// Make arealdistortion maps
process calc_areal_distortion {

    label 'connectome'
    stageInMode 'copy'

    publishDir "$params.out/sim_mesh/${sub}/registration/", \
                saveAs: { "${sub}.$it" }, \
                mode: 'copy'
    
    input:
    set val(sub), val(hemi), file(sphere), file(msm_sphere) from msm_outputs

    """
    wb_command -surface-distortion ${sphere} ${msm_sphere} \
                ${hemi}.areal_distortion.shape.gii
    """

}

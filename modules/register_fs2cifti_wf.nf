nextflow.preview.dsl=2


process convert_sulcal{
    /*
    Convert sulcal depth map to GIFTI

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere code
        sulc (Path): Sulcal depth file
        white (Path): Freesurfer white surface file

    Outputs:
        white_gifti (channel): (sub, white_gifti: Path) 
    */
        

    label 'freesurfer'
    containerOptions "-B ${params.license}:/license"

    input:
    tuple val(sub), val(hemi), path(sulc), path(white)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.sulc.native.shape.gii")

    """
    export FS_LICENSE=/license/license.txt
    mris_convert -c $sulc $white ${sub}.${hemi}.sulc.native.shape.gii
    """

}

process assign_sulcal{

    /*
    Assign structure to white sulcal depth GIFTI file

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        structure (str): Structure to assign to GIFTI file
        sulc (Path): Sulcal depth GIFTI file

    Output:
        assigned_sulc (channel): (sub, hemi, assigned_sulc: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), val(structure), path(sulc)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.assigned_sulc.native.shape.gii")

    """
    cp -L $sulc ${sub}.${hemi}.assigned_sulc.native.shape.gii
    wb_command -set-structure ${sub}.${hemi}.assigned_sulc.native.shape.gii $structure
    """


}

process invert_sulcal{
    /*
    Invert a sulcal depth map to a sulcal height map

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        sulc (Path): Sulcal depth GIFTI file

    Output:
        inverted_sulc (channel): (sub, hemi, inverted_sulc: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sulc)

    output:
    tuple val(sub), val(hemi), path("${sub}.inverted_sulc.${hemi}.shape.gii")

    """
    wb_command -metric-math 'a*(-1)' -var 'a' $sulc "${sub}.inverted_sulc.${hemi}.shape.gii"
    """
}

process convert_sphere{
    /*
    Convert a Freesurfer sphere to GIFTI format

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        sphere (Path): Input sphere to convert
        output (str): Name to assign to sphere

    Outputs:
        converted_sphere (channel): (sub, hemi, converted_sphere: Path)
    */

    label 'freesurfer'
    containerOptions "-B ${params.license}:/license"

    input:
    tuple val(sub), val(hemi), path(sphere), val(output)

    output:
    tuple val(sub), val(hemi), path("${sub}.${output}.surf.gii")

    shell:
    '''
    export FS_LICENSE=/license/license.txt
    mris_convert !{sphere} !{sphere}.surf.gii
    mv !{sphere}.surf.gii !{sub}.!{output}.surf.gii
    '''
}

process assign_sphere{

    /*
    Assign structure to registration sphere GIFTI file

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        structure (str): Structure to assign to GIFTI file
        sphere (Path): Registration sphere GIFTI file

    Output:
        assigned_sphere (channel): (sub, hemi, assigned_sphere: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), val(structure), path(sphere)

    output:
    tuple val(sub), val(hemi), path("${sub}.assigned_${sphere}")

    """
    cp -L ${sphere} ${sub}.assigned_${sphere}
    wb_command -set-structure ${sub}.assigned_${sphere} ${structure} -surface-type "SPHERICAL"
    """

}

process deform_sphere{

    /*
    Perform registration-sphere based deformation to 164k fsLR

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        sphere (path): Registation sphere to deform to 164k

    Outputs:
        reg_LR_sphere (channel): (sub, hemi, reg_LR: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sphere)

    output:
    tuple val(sub), val(hemi), file("${sub}.${hemi}.sphere.reg.reg_LR.native.surf.gii")

    shell:
    '''
    wb_command -surface-sphere-project-unproject \
                !{sphere} \
                /atlas/fsaverage.!{hemi}.sphere.164k_fs_!{hemi}.surf.gii \
                /atlas/fs_!{hemi}-to-fs_LR_fsaverage.!{hemi}_LR.spherical_std.164k_fs_!{hemi}.surf.gii \
                !{sub}.!{hemi}.sphere.reg.reg_LR.native.surf.gii
    '''

}

process spherical_affine{

    /*
    Perform affine alignment of registration spheres

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        sphere (Path): Registration sphere
        reg_LR_sphere (Path): Deformed registration sphere

    Outputs:
        affine (channel): (sub, hemi, affine: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sphere), path(reg_LR_sphere)

    output:
    tuple val(sub), val(hemi), path("${sub}_${hemi}_affine.mat")

    """
    wb_command -surface-affine-regression \
                ${sphere} \
                ${reg_LR_sphere} \
                ${sub}_${hemi}_affine.mat
    """
}

process normalize_rotation{
    /*
    Normalize rotation matrix to have det(A) = 1

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        affine (Path): Affine to normalize

    Outputs:
        normalized_affine (channel): (sub, hemi, norm_affine: Path)
    */
        

    label 'fieldopt'

    input:
    tuple val(sub), val(hemi), path(affine)

    output:
    tuple val(sub), val(hemi), path("${sub}_norm_affine.mat")

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
    np.savetxt("!{sub}_norm_affine.mat",M)
    '''

}

process apply_affine{
    /*
    Apply affine transformation to a sphere

    Arguments:
        sub (str): Subject ID
        hemi (Union['L' ,'R']): Hemisphere
        sphere (Path): Sphere to transform
        affine (Path): Normalized affine transformation

    Output:
        rot_sphere (channel): (sub, hemi, sphere_rot: Path)
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sphere), path(affine)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.sphere_rot.surf.gii")

    """
    wb_command -surface-apply-affine \
                $sphere \
                $affine \
                ${sub}.${hemi}.sphere_rot.surf.gii

    wb_command -surface-modify-sphere \
                ${sub}.${hemi}.sphere_rot.surf.gii \
                100 \
                ${sub}.${hemi}.sphere_rot.surf.gii
    """

}

process msm_sulc{
    /*
    Perform MSM-Sulcal depth-based registration to 164k fsLR sphere

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere
        sphere (Path): Registration sphere 
        sulc (Path): Sulcal height map
        structure (str): Structure to assign to resultant sphere

    Outputs:
        reg_msm (channel): (sub, hemi, sphere, reg_msm: Path) MSM registration sphere 
    */


    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sphere), path(sulc), val(structure)

    output:
    tuple val(sub), val(hemi), path(sphere), path("${sub}.${hemi}.sphere.reg_msm.surf.gii")

    shell:
    '''
    /msm/msm --inmesh=!{sphere} \
             --indata=!{sulc} \
             --refmesh=/atlas/fsaverage.!{hemi}_LR.spherical_std.164k_fs_LR.surf.gii \
             --refdata=/atlas/!{hemi}.refsulc.164k_fs_LR.shape.gii \
             --conf=/msm_conf/MSMSulcStrainFinalconf \
             --out=!{hemi}. \
             --verbose

    mv "!{hemi}.sphere.reg.surf.gii" \
       "!{sub}.!{hemi}.sphere.reg_msm.surf.gii"

    wb_command -set-structure !{sub}.!{hemi}.sphere.reg_msm.surf.gii \
                                !{structure}
    '''

}

process areal_distortion{

    /*
    Compute areal surface distortion map from an MSM registration

    Arguments:
        sub (str): Subject ID        
        hemi (Union['L', 'R']): Hemisphere
        sphere (Path): Input registration sphere
        msm_sphere (Path): MSM registration sphere

    Outputs:
        areal_distortion (channel): (sub, hemi, areal: Path)
    */
        

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(sphere), path(msm_sphere)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.areal_distortion.shape.gii"), emit: areal

    """
    wb_command -surface-distortion \
                ${sphere} \
                ${msm_sphere} \
                ${sub}.${hemi}.areal_distortion.shape.gii
    """




}

workflow registration_wf {
    /*
    Perform MSM surface-based registration between Freesurfer native sphere to fs_LR 164k

    Arguments:
        fs_dirs (channel): (sub, fs_dir: Path) Freesurfer directory for subject

    Outputs:
        msm_sphere (channel): (sub, hemi: Union['L', 'R'], msm_sphere: Path) Hemispheric MSM registration sphere
    */


    take:
        fs_dirs

    main:

        // Might have to migrate this over fs2gifti
        // Convert sulcal information from freesurfer to connectome workbench
        sulcal_input = fs_dirs
                            .spread ( ['L','R'] )
                            .map{s,f,h ->   [
                                                s,
                                                h,
                                                "${f}/surf/${h.toLowerCase()}h.sulc",
                                                "${f}/surf/${h.toLowerCase()}h.white"
                                            ]
                                }
        convert_sulcal(sulcal_input)

        // Assign structure to sulcal map then invert
        structure_map = ['L' : 'CORTEX_LEFT', 'R' : 'CORTEX_RIGHT' ]
        assign_input = convert_sulcal.out
                                    .map{ s,h,g ->  [
                                                        s,h,
                                                        structure_map[h],
                                                        g
                                                    ]
                                        }
        assign_sulcal(assign_input)
        invert_sulcal(assign_sulcal.out)

        // Now convert spheres over, assign properties,
        registration_spheres = fs_dirs
                                    .spread( ['L','R'] )
                                    .spread( ['sphere','sphere.reg'] )
                                    .map{ s,f,h,sph ->  [
                                                            s,h,
                                                            "${f}/surf/${h.toLowerCase()}h.${sph}",
                                                            "${h}.${sph}"
                                                        ]
                                        }
        convert_sphere(registration_spheres)
        assign_sphere_input = convert_sphere.out
                                        .map{ s,h,sph ->[
                                                            s,h,
                                                            structure_map[h],
                                                            sph
                                                        ]
                                            }
        assign_sphere(assign_sphere_input)

        //Pull reg spheres and perform spherical deformation
        deform_sphere_input = assign_sphere.out
                                    .filter { it[2].name.contains('reg') }
        deform_sphere(deform_sphere_input)

        // Merge with native sphere and compute affine
        affine_input = assign_sphere.out
                                    .filter { !(it[2].name.contains('reg')) }
                                    .join(deform_sphere.out, by : [0,1])
        spherical_affine(affine_input)

        // Normalize affine transformation
        normalize_rotation(spherical_affine.out)

        // Apply affine transformation to sphere
        rotation_input = assign_sphere.out
                                .join(normalize_rotation.out, by: [0,1])

        apply_affine(rotation_input)


        // Perform MSM
        msm_input = apply_affine.out
                            .join(invert_sulcal.out, by: [0,1])
                            .map{ s,h,sph,sulc ->   [
                                                        s,h,sph,sulc,
                                                        structure_map[h]
                                                    ]
                                }
        msm_sulc(msm_input)

        // Make areal distortion map
        areal_distortion(msm_sulc.out)

        // Make outputs
        msm_sphere_out = msm_sulc.out
                                .map{ s,h,sph,msm -> [s,h,msm] }

        emit:
            msm_sphere = msm_sphere_out

}

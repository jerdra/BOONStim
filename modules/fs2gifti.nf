nextflow.preview.dsl=2


process convert_fs2gifti{

    /*
    Convert Freesurfer surface files to GIFTI format

    Arguments:
        sub (str): Subject ID
        hemi (Union['L','R']): Hemisphere key
        surf (str): Freesurfer surface type
        fs_surf (Path): Path to surface

    Parameters:
        license (Path): Freesurfer license file

    Outputs:
        hemi_surf (channel): (subject, hemi: Union['L', 'R'], surf: Path) GIFTI surf file
    */

    label 'freesurfer'
    maxForks 4

    containerOptions "-B ${params.license}:/license"

    input:
    tuple val(sub), val(hemi), val(surf), path(fs_surf)

    output:
    tuple val(sub), val(hemi), val(surf),\
    path("${sub}.${hemi}.unassigned_${surf}.surf.gii"), emit: hemi_surf

    shell:
    '''
    export FS_LICENSE=/license/license.txt

    mris_convert !{fs_surf} !{surf}.surf.gii

    hemi="!{hemi}"
    mv ${hemi,,}h.!{surf}.surf.gii \
        !{sub}.!{hemi}.unassigned_!{surf}.surf.gii
    '''

}


process assign_structure {

    /*
    Assign workbench structure to GIFTI file

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere key
        surf (str): Type of surface
        structure (str): Surface type to assign to GIFTI
        gifti (Path): Path to GIFTI file

    Outputs:
        gifti (channel): (sub, hemi: Union['L','R'], gifti: Path) Path to
            GIFTI file with assigned `structure`
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), val(surf), val(structure), path(gifti)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.${surf}.surf.gii"), emit: gifti

    shell:
    '''
    cp -L !{gifti} !{sub}.!{hemi}.!{surf}.surf.gii
    wb_command -set-structure !{sub}.!{hemi}.!{surf}.surf.gii !{structure}

    '''

}

process compute_midthickness {
    /*
    Calculate the midthickness surface of pial/white surface

    Arguments:
        sub (str): Subject ID
        hemi (Union['L', 'R']): Hemisphere key
        pial (Path): Path to pial surface
        white (Path): Path to white surface

    Outputs:
        midthickness (Path): Path to computed midthickness file
    */

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(pial), path(white)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.midthickness.surf.gii"), emit: midthickness

    """
    wb_command -surface-average -surf ${pial} \
                                -surf ${white} \
                                ${sub}.${hemi}.midthickness.surf.gii
    """

}


workflow make_giftis {

    /*
    Convert Freesurfer pial, white surface files into GIFTI format. A
    midthickness GIFTI file is also calculated

    Arguments:
        fs_dirs (channel): (subject, fs_dir: Path) Path to subject Freesurfer directory

    Parameters:
        license (Path): Path to Freesurfer license file

    Outputs:
        midthickness (channel): (subject, hemi: Union['L', 'R'], midthick: Path)
            Midthickness GIFTI file
        pial (channel): (subject, hemi: Union['L', 'R'],  pial: Path) Pial GIFTI files
        white (channel): (subject, hemi: Union['L', 'R'],  white: Path) White GIFTI files
    */



    take:
    fs_dirs

    main:
        // Construct pathing
        sub_surfaces = fs_dirs
                            .spread(['L','R'])
                            .spread(['pial','white'])
                            .map{ n,fs,h,surf ->   [
                                                        n,h,surf,
                                                        "${fs}/surf/${h.toLowerCase()}h.${surf}"
                                                   ]
                                }
        convert_fs2gifti(sub_surfaces)

        // Map hemisphere over to structure
        structure_map = ['L' : 'CORTEX_LEFT', 'R' : 'CORTEX_RIGHT' ]
        assign_structure_input = convert_fs2gifti.out
                                                .map{ s,h,surf,g ->  [
                                                                    s,
                                                                    h,
                                                                    surf,
                                                                    structure_map[h],
                                                                    g
                                                                ]
                                                    }
        assign_structure(assign_structure_input)

        // Group together surfaces to subject/hemisphere
        midthickness_input = assign_structure.out
                                        .groupTuple(by: [0,1], sort: {it.baseName}, size: 2)
                                        .map{ n,h,surfs ->  [
                                                                n,h,
                                                                surfs[0],surfs[1]
                                                            ]
                                            }
        compute_midthickness(midthickness_input)

        pial_out = midthickness_input.map { n,h,p,w -> [n,h,p] }
        white_out = midthickness_input.map { n,h,p,w -> [n,h,w] }

        emit:
            midthickness = compute_midthickness.out.midthickness
            pial = pial_out
            white = white_out


}

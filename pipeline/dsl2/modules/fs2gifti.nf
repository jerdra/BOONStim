nextflow.preview.dsl=2


process convert_fs2gifti{

    label 'freesurfer'
    maxForks 4

    containerOptions "-B ${params.license}:/license"

    input:
    tuple val(sub), val(hemi), val(surf), path(fs_surf)

    output:
    tuple val(sub), val(hemi), path("${hemi}.${surf}.surf.gii"), emit: hemi_surf

    shell:
    '''
    export FS_LICENSE=/license/license.txt

    mris_convert !{fs_surf} !{hemi}.!{surf}.surf.gii
    
    hemi="!{hemi}"
    mv ${hemi,,}h.!{hemi}.!{surf}.surf.gii \
        !{hemi}.!{surf}.surf.gii
    '''

}


process assign_structure {

    label 'connectome'

    input:
    tuple val(sub), val(hemi), val(structure), path(gifti)

    output:
    tuple val(sub), val(hemi), path("assigned_$gifti"), emit: gifti

    shell:
    '''
    cp -L !{gifti} assigned_!{gifti}
    wb_command -set-structure assigned_!{gifti} !{structure}
    '''

}

process compute_midthickness {

    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(pial), path(white)

    output:
    tuple val(sub), val(hemi), path("${hemi}.midthickness.surf.gii"), emit: midthickness

    """
    wb_command -surface-average -surf ${pial} \
                                -surf ${white} \
                                ${hemi}.midthickness.surf.gii
    """

}


workflow make_giftis {
    
    get: fs_dirs

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
                                                .map{ s,h,g ->  [
                                                                    s,
                                                                    h,
                                                                    structure_map[h],
                                                                    g
                                                                ]
                                                    }
        assign_structure(assign_structure_input)

        // Group together surfaces to subject/hemisphere
        midthickness_input = assign_structure.out
                                        .groupTuple(by: [0,1], sort: {it.baseName})
                                        .map{ n,h,surfs ->  [
                                                                n,h,
                                                                surfs[0],surfs[1]
                                                            ]
                                            }
        compute_midthickness(midthickness_input)

        emit:
            midthickness = compute_midthickness.out.midthickness
            surfaces = midthickness_input
            
            
}

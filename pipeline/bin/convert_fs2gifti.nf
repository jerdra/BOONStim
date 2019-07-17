
// Nextflow pipeline to convert Freesurfer files to GIFTI format

//INPUTS:
//out               Output directory containing mri2mesh outputs



//OUTPUTS:
//Script will output to output_dir/sim_mesh/sub/surfaces


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


//Get all surfaces
sub_surfaces = sub_channel
                    .map { n -> [n.baseName, n] }
                    .spread(['L','R'])
                    .spread(['pial','white'])
                    .map { n -> [
                                    n[0],
                                    n[2],
                                    n[3],
                                    n[1] + "/fs_${n[0]}/surf/${n[2].toLowerCase()}h.${n[3]}"
                                ]
                        }


structure_map = ['L': 'CORTEX_LEFT',
                 'R': 'CORTEX_RIGHT']

// Freesurfer --> Gifti
process convert_fs2gifti {

    label 'freesurfer'
    maxForks 4
    stageInMode "copy"

    containerOptions "-B ${params.license}:/license" 

    input:
    set val(sub), val(hemi), val(surf), file(fs_surf) from sub_surfaces

    output:
    set val(sub), val(hemi), file("${hemi}.${surf}.surf.gii") into sub_gifti

    shell:
    '''
    set +u
    export FS_LICENSE=/license/license.txt
    
    mris_convert !{fs_surf} !{hemi}.!{surf}.surf.gii

    #Cause freesurfer wants to prepend fs hemisphere
    hemi="!{hemi}"
    mv ${hemi,,}h.!{hemi}.!{surf}.surf.gii \
        !{hemi}.!{surf}.surf.gii
    '''
    
}

assign_struct_input = sub_gifti
                            .map { n -> [
                                            n[0],
                                            n[1],
                                            structure_map[n[1]],
                                            n[2]
                                        ]
                                }


// Assign connectome workbench structures
process assign_structure {

    label 'connectome'

    input:
    set val(sub), val(hemi), val(structure), file(gifti) from assign_struct_input

    output:
    set val(sub), val(hemi), file(gifti) into assigned_gifti

    """
    wb_command -set-structure ${gifti} ${structure}
    """

}



//OUTPUT: [sub, hemi, pial, white]
midthickness_input = assigned_gifti
                            .groupTuple(by: [0,1],
                                        sort: { it.baseName } 
                                    )
                            .map { n -> [
                                            n[0],
                                            n[1],
                                            n[2][0],
                                            n[2][1]
                                        ]
                                }

//Compute midthickness
process compute_midthickness {

    label 'connectome'

    input:
    set val(sub), val(hemi), file(pial), file(white) from midthickness_input

    output:
    set val(sub), val(hemi), file("${hemi}.midthickness.surf.gii") into midthickness

    """
    wb_command -surface-average -surf ${pial} \
                -surf ${white} \
                ${hemi}.midthickness.surf.gii
    """


}
midthickness.subscribe{ log.info("$it") }

nextflow.preview.dsl = 2

process adm_optimization {
/*
*   Inputs:
*   sub - Subject name
*   msh - Path to gmsh .msh file
*   target_coord - Path to numpy file containing (X,Y,Z) coordinates to target
*   target_radius - ROI radius
*   coil_distance - Distance of coil from head
*   coil - Path to coil definition file
*/

    input:
    tuple val(sub), path(msh),\
    path(target_coord), val(target_radius),\
    val(coil_distance), path(coil)

    output:
    tuple val(sub), path("${sub}_optimized_fields.msh"), emit: fields
    tuple val(sub), path("${sub}_optimized_coil.geo"), emit: coil
    tuple val(sub), path("${sub}_orientation.npy"), emit: coords

    shell:
    '''
    /bin/adm.py !{msh} \
    !{target_coord} !{target_radius} \
    !{coil_distance} !{coil} \
    !{sub}_orientation.npy \
    --msh !{sub}_optimized_fields.msh \
    --geo !{sub}_optimized_geo.geo
        
    '''
}

workflow optimize_wf {
    take:
        msh
        target_coord
        target_radius
        coil_distance
        coil

    main:
        adm_optimization(
            msh
                .join(target_coord)
                .spread([target_radius])
                .spread([coil_distance])
                .spread([coil])
        )
                
    emit:
        coil = adm_optimization.out.coil
        fields = adm_optimization.out.fields
        coords = adm_optimization.out.coords

        
}

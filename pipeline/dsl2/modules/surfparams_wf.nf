nextflow.preview.dsl=2

process extract_surf_patch {

    label 'numpy'
    containerOptions "-B ${params.bin}:/scripts"
    
    input:
    tuple val(sub), path(msh), path(centroid)
    
    output:
    tuple val(sub), path('*dilated_coords.npy'), path('*mean_norm.npy'), emit: surf_patch

    shell:
    '''
    /scripts/extract_surface_patch.py !{msh} !{centroid} "surf"
    '''
}

process parameterize_surf {

    label 'numpy'
    containerOptions "-B ${params.bin}:/scripts"
    
    input:
    tuple val(sub), path(patch), path(norm)

    output:
    tuple val(sub), path('*_C.npy'), emit: C
    tuple val(sub), path('*_R.npy'), emit: R
    tuple val(sub), path('*_bounds.npy'), emit: bounds
    
    shell:
    '''
    /scripts/parameterize_surface_patch.py !{patch} !{norm} "out"
    '''
}


workflow parameterization_wf {

    get:
        msh 
        centroid

    main:
        
        // Make surface patch for fitting
        surf_patch_input = msh.join(centroid, by: 0)
        extract_surf_patch(surf_patch_input)
        extract_surf_patch.out.surf_patch

        // Parameterize surface patch
        parameterize_surf(extract_surf_patch.out.surf_patch)

    emit:
        C = parameterize_surf.out.C
        R = parameterize_surf.out.R
        bounds = parameterize_surf.out.bounds

}

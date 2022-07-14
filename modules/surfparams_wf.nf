nextflow.preview.dsl=2

process extract_surf_patch {

    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(msh), path(centroid)

    output:
    tuple val(sub), path("${sub}_dilated_coords.npy"), path("${sub}_mean_norm.npy"), emit: surf_patch
    tuple val(sub), path("${sub}_param_surf.msh"), emit: qc_surf

    shell:
    '''
    /scripts/extract_surface_patch.py !{msh} !{centroid} !{sub}
    '''
}

process parameterize_surf {

    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(patch), path(norm)

    output:
    tuple val(sub), path("${sub}_C.npy"), emit: C
    tuple val(sub), path("${sub}_R.npy"), emit: R
    tuple val(sub), path("${sub}_bounds.npy"), emit: bounds

    shell:
    '''
    /scripts/parameterize_surface_patch.py !{patch} !{norm} !{sub}
    '''
}

process qc_parameterization {

    label 'fieldopt'
    containerOptions "-B ${params.bin}:/scripts"

    input:
    tuple val(sub), path(qc_surf), path(C), path(R), path(bounds)

    output:
    tuple val(sub), path("${sub}_quadratic_surf.msh"), emit: qc_param

    shell:
    '''
    /scripts/qc_parametric.py !{qc_surf} !{C} !{R} !{bounds} !{sub}_quadratic_surf.msh
    '''
}


workflow parameterization_wf {

    take:
        msh
        centroid

    main:

        // Make surface patch for fitting
        surf_patch_input = msh.join(centroid, by: 0)
        extract_surf_patch(surf_patch_input)

        // Parameterize surface patch
        parameterize_surf(extract_surf_patch.out.surf_patch)

        qc_input = extract_surf_patch.out.qc_surf
                                     .join(parameterize_surf.out.C, by: 0)
                                     .join(parameterize_surf.out.R, by: 0)
                                     .join(parameterize_surf.out.bounds, by: 0)
        qc_parameterization(qc_input)

    emit:
        C = parameterize_surf.out.C
        R = parameterize_surf.out.R
        bounds = parameterize_surf.out.bounds
        qc_param = qc_parameterization.out.qc_param

}

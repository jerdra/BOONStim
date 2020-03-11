nextflow.preview.dsl=2

include calculate_weightfunc_wf from './compute_dmpfc_connectivity.nf' params(params)
include mask_wf from './make_dmpfc_mask.nf' params(params)

process threshold_weightfunc{

    label "connectome"

    input:
    tuple val(sub), path(wf), path(mask)

    output:
    tuple val(sub), path("${sub}.thresholded_weightfunc.dscalar.nii"), emit: weightfunc

    shell:
    '''

    # Get 80th percentile in mask
    thres=$(wb_command -cifti-stats "masked.dscalar.nii" -percentile 80 \
            -roi !{mask})
    wb_command -cifti-math \
                "x * $thres" \
                -var "x" !{wf} \
                "!{sub}.thresholded_weightfunc.dscalar.nii"
    '''


}

workflow weightfunc_wf {

    get: derivatives

    main:
        cifti = derivatives.map{s,f,c -> [s,c]}
        calculate_weightfunc_wf(derivatives)
        mask_wf(cifti)

        threshold_input = calculate_weightfunc_wf.out.weightfunc
                                    .join(mask_wf.out.mask)
        threshold_weightfunc(threshold_input)

    emit:
        mask = mask_wf.out.mask
        weightfunc = threshold_weightfunc.out.weightfunc
}

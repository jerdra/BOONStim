nextflow.preview.dsl=2

include calculate_weightfunc_wf from './compute_dmpfc_connectivity.nf' params(params)
include mask_wf from './make_dmpfc_mask.nf' params(params)

workflow weightfunc_wf {

    get: derivatives

    main:

        cifti = derivatives.map { s,f,c -> [s,c] }
        calculate_weightfunc_wf(derivatives)
        mask_wf(cifti,calculate_weightfunc_wf.out.weightfunc)

    emit:
        mask = mask_wf.out.mask
        weightfunc = calculate_weightfunc_wf.out.weightfunc
        
        



}

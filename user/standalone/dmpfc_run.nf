nextflow.preview.dsl=2

if (!params.derivatives){

    log.info('Insufficient input specification!')
    log.info('Needs --derivatives!')
    System.exit(1)
}

if (params.subjects) {

    log.info("Subject list file provided: $params.subjects")

}

// Import weightfunc workflow
include weightfunc_wf from '../calculate_mentalizing_weightfunc.nf' params(params)

if (params.subjects){
    sublist = file(params.subjects)
    bids_channel = Channel.from(sublist)
                            .splitTest() { it.strip() }
}else{
    ciftify_dir="$params.derivatives/ciftify"
    all_subs = file(ciftify_dir).list()
    bids_channel = Channel.from(all_subs)
                            .filter { it.contains('sub-') }
}

workflow {

    main:

        // Using derivatives arguments, get...
        ciftify_dirs = Channel.fromPath("$params.derivatives/ciftify/sub-*", type: 'dir')
                            .map{ c -> [c.getFileName().getName(),c]}
        fmriprep_dirs = Channel.fromPath("$params.derivatives/fmriprep/sub-*", type: 'dir')
                            .map{ f -> [f.getFileName().getName(),f]}
        input_channel = bids_channel.join(fmriprep_dirs).join(ciftify_dirs)

        // Feed into pipeline
        weightfunc_wf(input_channel)

    publish:
        weightfunc_wf.out.mask to: "$params.out", mode: 'copy'
        weightfunc_wf.out.weightfunc to: "$params.out", mode: 'copy'



}

// BOONSTIM FULL PIPELINE WORKFLOW
nextflow.preview.dsl=2

if (!params.bids || !params.out){

    log.info('Insufficient input specification!')
    log.info('Needs --bids, --out!')
    log.info('Exiting...')
    System.exit(1)

}

log.info("BIDS Directory: $params.bids")
log.info("Output Directory: $params.out")

//Subject list flag
if (params.subjects) {
    log.info ("Subject file provided: $params.subjects")
}

//// Extract subjects to run
all_dirs = file(params.bids).list()
input_dirs = new File(params.bids).list()
output_dirs = new File(params.out).list()

if (params.subjects) {
    sublist = file(params.subjects)
    bids_channel = Channel.from(sublist)
                               .splitText() { it.strip() }
} else{
    bids_channel = Channel.from(all_dirs)
}


process run_boonstim_subject{

    module 'slurm'

    input:
    val(sub)
    
    output:
    val(sub)

    shell:
    '''
    #!/bin/bash
    
    nextflow !{params.boonstim}/boonstim.nf \
                -c !{params.boonstim}/config/boonstim.nf.config \
                --bids !{params.bids} \
                --out !{params.out} \
                --subject !{sub} \
                -profile local
    '''
}

workflow{
    run_boonstim_subject(bids_channel)
}

/*

Run optimization routine given a derivatives directory containing full input
specification for each subject.

Uses fem_optimization derivatives directory

*/


if (!params.out){

    log.info('Insufficient input specification!')
    log.info('Needs --out!')
    log.info('Exiting...')
    System.exit(1)

}


// Get a list of subjects with the required output
fem_subs = "$params.out/sim_mesh/sub-*/fem_optimization"

//Get full input requirements
input_subs = Channel.fromPath(fem_subs, type: 'dir')
                    //.subscribe{ log.info("$it") }
                    .map    { n ->  [
                                        n.getParent().getBaseName(),
                                        n
                                    ]
                            }
                    .map    { n,p-> [
                                        n,
                                        file("$p/${n}_tetraweights.npy"),
                                        file("$p/${n}_surf_C.npy"),
                                        file("$p/${n}_surf_R.npy"),
                                        file("$p/${n}_surf_bounds.npy"),
                                        file("$params.out/sim_mesh/${n}/${n}.msh")
                                    ]
                            }
                    .subscribe { log.info("$it") }

//Run optimization routine
process fem_optimize{

    input:
    set val(sub),
        file(W),
        file(C), file(R), file(b),
        file(msh) \
    from input_subs
    

    shell:
    '''
    !{params.rtms_bin}/optimize_fem.py !{msh} \
                        !{W} !{C} !{b} !{R} \
                        !{params.coil} \
                        loc.txt rot.txt \
                        --cpus !{params.cpus}
    '''

}


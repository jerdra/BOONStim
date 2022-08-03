nextflow.preview.dsl = 2

process get_scalp_seed {

    /*
    Compute a good candidate seed position given a scalar-valued map

    Arguments:
        sub (str): Subject ID
        mesh (Path): Realistic head model
        dscalar (Path): Scalar-valued map
        l_pial (Path): Left pial GIFTI surface
        r_pial (Path): Right pial GIFTI surface

    Outputs:
        seed (channel): (sub, seed: Path) Candidate coordinate on subject scalp
        qchtml (channel): (sub, html: Path) QC page
    */

    label 'fieldopt'

    input:
    tuple val(sub), path(mesh), path(dscalar), path(l_pial), path(r_pial)

    output:
    tuple val(sub), path("${sub}_seed.txt"), emit: seed
    tuple val(sub), path("${sub}_qcseed.html"), emit: qchtml

    shell:
    '''
    #!/bin/bash

    /scripts/get_scalp_seed.py !{mesh} !{dscalar} !{l_pial} !{r_pial} \
                               !{sub}_seed.txt --qc-file !{sub}_qcseed.html
    '''
}

def lr_branch = branchCriteria {
                left: it[1] == 'L'
                    return [it[0], it[2]]
                right: it[1] == 'R'
                    return [it[0],it[2]]
                }

workflow centroid_radial_wf{

    /*
    Compute approximate radial projection out from pial surface onto scalp
    at a given coordinate

    Arguments:
        msh (channel): (subject, msh: Path) Subject .msh file
        dscalar (channel): (subject, dscalar: Path) dscalar file to compute projection from
        pial (channel): (subject, hemi: Union['L' 'R'], pial: Path) Pial file

    Outputs:
        centroid (channel): (subject, centroid: Path) Computed centroid
        qc (channel): (subject, qchtml: Path) QC file associated with `centroid`
    */


    take:
        msh
        dscalar
        pial

    main:

        // Branch out left/right surfaces
        pial.branch(lr_branch).set{pial_surfs}

        i_get_scalp_seed = msh.join(dscalar)
                              .join(pial_surfs.left)
                              .join(pial_surfs.right)

        get_scalp_seed(i_get_scalp_seed)

    emit:
        centroid = get_scalp_seed.out.seed
        qc = get_scalp_seed.out.qchtml


}

nextflow.preview.dsl=2

include { weightfunc_wf } from "${params.weightworkflow}" params(params)
include { cifti_meshing_wf as cifti_mesh_wf } from '../modules/cifti_mesh_wf.nf' params(params)
include { rigidRegistration, coordinate_transform, map_coordinate } from "../modules/transformation.nf"

// Default params
params.radius = 20 // in mm


process format_for_ants {

    /*
    Format a RAS coordinate into an LPS ANTS CSV file

    Arguments:
        sub (str): Subject ID
        coords (Path): Path to coordinate .txt file

    Output:
        ants_csv (channel): (sub, csv: Path) Ants formatted CSV file
    */

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_pretransform.txt"), emit: ants_csv

    shell:
    """
    echo "x,y,z,t" > "${sub}_pretransform.txt"
    echo "${x},${y},${z},0" >> "${sub}_pretransform.txt"
    """
}

process format_from_ants {

    /*
    Format an  LPS ANTS CSV file into a RAS coordinate file

    Arguments:
        sub (str): Subject ID
        coords (Path): Path to ants coordinate .txt file

    Output:
        ras_coords (channel): (sub, npy: Path) Path to .npy file in RAS
    */

    input:
    tuple val(sub), path(coords)

    output:
    tuple val(sub), path("${sub}_ras_coords.npy"), emit: ras_coords

    shell:
    """
    #!/usr/bin/env python

    import numpy as np

    # LPS space
    coords = np.loadtxt('${sub}_ras_coords.txt')
    coords[0] = -coords[0]
    coords[1] = -coords[1]
    np.savetxt("${sub}_ras_coords.txt", coords[:-1], delim=',')
    """
}

workflow coordinate_optimization {

    /*
    Perform single-coordinate optimization with a set radius

    Arguments:
        subject_channel (channel): Subject IDs

    Parameters:
        radius (float): Radius of optimization target around derived coordinate
    */

    take:
    subject_channel

    main:
    cifti_meshing_wf(subject_channel)

    weightfunc_input = cifti_mesh_wf.out.fmriprep
                                        .join( cifti_mesh_wf.out.cifti, by : 0 )
    weightfunc_input
    weightfunc_wf(weightfunc_input)

    i_rigidRegistration = cifti_meshing_wf.out.cifti
        .join(cifti_meshing_wf.out.mesh_m2m, by: 0)
        .map{ s, c, m2m -> [
            s,
            "${c}/T1w/T1w.nii.gz",
            "${m2m}/T1fs_resamp.nii.gz"
            ]
        }


    rigidRegistration(i_rigidRegistration)
    map_coordinate(weightfunc_wf.out.coordinate, rigidRegistration.out.transform)

    i_adm = cifti_meshing_wf.out.msh
        .join(cifti_meshing_wf.out.mesh_fs)
        .join(map_coordinate.out.ras_coords)
        .spread([params.radius])
    adm_wf(i_adm)
}

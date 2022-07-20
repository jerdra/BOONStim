nextflow.preview.dsl=2

include {weightfunc_wf} from "${params.weightworkflow}" params(params)
include {cifti_meshing_wf as cifti_mesh_wf} from '../modules/cifti_mesh_wf.nf' params(params)

workflow coordinate_optimization {

    take:
    subject_channel

    main:
    cifti_meshing_wf(subject_channel)

    weightfunc_input = cifti_mesh_wf.out.fmriprep
                                        .join( cifti_mesh_wf.out.cifti, by : 0 )
    weightfunc_input
    weightfunc_wf(weightfunc_input)



}

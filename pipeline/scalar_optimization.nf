nextflow.preview.dsl=2

include {weightfunc_wf} from "${params.weightworkflow}" params(params)
include {cifti_meshing_wf as cifti_mesh_wf} from '../modules/cifti_mesh_wf.nf' params(params)
include {make_giftis} from '../modules/fs2gifti.nf' params(params)
include {registration_wf} from '../modules/register_fs2cifti_wf.nf' params(params)
include {resample2native_wf as resamplemask_wf} from '../modules/resample2native.nf' params(params)
include {resample2native_wf as resampleweightfunc_wf} from '../modules/resample2native.nf' params(params)
include {resample2native_wf as resampledistmap_wf} from '../modules/resample2native.nf' params(params)
include {resample2native_wf as resamplesulc_wf} from '../modules/resample2native.nf' params(params)
include {centroid_radial_wf as centroid_wf} from '../modules/centroid_wf.nf' params(params)
include {tet_project_wf as tet_project_weightfunc_wf} from '../modules/tetrahedral_wf.nf' params(params)
include {tet_project_wf as tet_project_roi_wf} from '../modules/tetrahedral_wf.nf' params(params)
include {calculate_reference_field_wf} from '../modules/reference_field_wf.nf' params(params)
include {fieldscaling_wf} from '../modules/field_scaling.nf' params(params)
include {optimize_wf} from "../modules/optimization.nf" params(params)
include {apply_mask as centroid_mask} from '../modules/utils.nf' params(params)
include {apply_mask as weightfunc_mask} from '../modules/utils.nf' params(params)
include {cifti_dilate as dilate_mask} from '../modules/utils.nf' params(params)

workflow scalar_optimization {
    take:
        subject_channel

    main:

        // Main preprocessing routine
        cifti_mesh_wf(subject_channel)

        make_giftis(cifti_mesh_wf.out.mesh_fs)
        registration_wf(cifti_mesh_wf.out.mesh_fs)

        // User-defined weightfunction workflow
        weightfunc_input = cifti_mesh_wf.out.fmriprep
                                            .join( cifti_mesh_wf.out.cifti, by : 0 )
        weightfunc_input
        weightfunc_wf(weightfunc_input)

        // Calculate centroid on resampled data
        centroid_mask_input = weightfunc_wf.out.weightfunc
                                            .join(weightfunc_wf.out.mask)
        centroid_mask(centroid_mask_input)
        resamplemask_wf(centroid_mask.out.masked, registration_wf.out.msm_sphere)

        // Tetrahedral workflow

        // Resample the weightfunction
        dilate_mask_input = weightfunc_wf.out.mask
                                            .join(cifti_mesh_wf.out.cifti, by: 0)
                                            .map{ s,w,c ->  [
                                                                s,w,
                                                                "${c}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                                                "${c}/MNINonLinear/fsaverage_LR32k/${s}.R.midthickness.32k_fs_LR.surf.gii"
                                                            ]
                                                }
        dilate_mask(dilate_mask_input)
        weightfunc_mask_input = weightfunc_wf.out.weightfunc
                                            .join(dilate_mask.out.dilated, by: 0)
        weightfunc_mask(weightfunc_mask_input)
        resampleweightfunc_wf(weightfunc_mask.out.masked, registration_wf.out.msm_sphere)

        // Calculate a scalp seed
        centroid_wf(cifti_mesh_wf.out.msh,
                    resampleweightfunc_wf.out.resampled,
                    make_giftis.out.pial)

        tet_project_weightfunc_wf(resampleweightfunc_wf.out.resampled,
                        make_giftis.out.pial,
                        make_giftis.out.white,
                        make_giftis.out.midthickness,
                        cifti_mesh_wf.out.t1fs_conform,
                        cifti_mesh_wf.out.msh)

        // Gather inputs for optimization
        optimize_wf(
                    cifti_mesh_wf.out.msh,
                    tet_project_weightfunc_wf.out.fem_weights,
                    centroid_wf.out.centroid,
                    params.coil
                   )

        // Calculate scaling factor between coil and cortex across multiple references
        calculate_reference_field_wf(cifti_mesh_wf.out.cifti,
                                     Channel.from(params.ref_coords))
        resampledistmap_wf(
            calculate_reference_field_wf.out.rois,
            registration_wf.out.msm_sphere
        )

        fieldscaling_wf(
                        optimize_wf.out.fields,
                        make_giftis.out.pial,
                        resampledistmap_wf.out.resampled,
                        optimize_wf.out.matsimnibs
                      )

        // Gather BOONStim outputs for publishing
        registration_wf.out.msm_sphere.branch(lr_branch).set { msm }
        make_giftis.out.pial.branch(lr_branch).set { pial }
        make_giftis.out.white.branch(lr_branch).set { white }
        make_giftis.out.midthickness.branch(lr_branch).set { midthick }

        /* Step 1: Publish base outputs */
        i_publish_base = cifti_mesh_wf.out.t1fs_conform
                                .join(centroid_wf.out.centroid)
                                .join(tet_project_weightfunc_wf.out.fem_weights)
                                .join(resampleweightfunc_wf.out.resampled)
        publish_base(i_publish_base)

        /* Step 2: Publish native space surfaces used to map out
        weight function
        */
        i_publish_surfs = pial.left.join(pial.right)
                            .join(white.left).join(white.right)
                            .join(midthick.left).join(midthick.right)
                            .join(msm.left).join(msm.right)
        publish_surfs(i_publish_surfs)

        /* Step 3: Publish meshing results from mri2mesh */
        i_publish_mri2mesh = cifti_mesh_wf.out.mesh_m2m
                                .join(cifti_mesh_wf.out.mesh_fs)
        publish_mri2mesh(i_publish_mri2mesh)

        /* Step 4: Publish optimization results */
        i_publish_opt = optimize_wf.out.fields
                            .join(optimize_wf.out.coil)
                            .join(optimize_wf.out.history)
                            .join(optimize_wf.out.brainsight)
                            .join(optimize_wf.out.localite)
        publish_opt(i_publish_opt)

        /* Step 5: Publish the reference scaling values */
        i_publish_scaleref = fieldscaling_wf.out.scaling_factor
                                            .join(fieldscaling_wf.out.qc_html, by: [0,1])
                                            .join(fieldscaling_wf.out.qc_geo, by: [0,1])
        publish_scaleref(i_publish_scaleref)

        // Publish Ciftify outputs
        publish_cifti_input = cifti_mesh_wf.out.cifti
                                    .join(cifti_mesh_wf.out.cifti_qc_fmri)
                                    .join(cifti_mesh_wf.out.cifti_qc_recon)
                                    .join(cifti_mesh_wf.out.fmriprep)
                                    .join(cifti_mesh_wf.out.fmriprep_html)
                                    .join(cifti_mesh_wf.out.freesurfer)
                                    .combine(["$params.zz"])
        publish_cifti(publish_cifti_input)
}
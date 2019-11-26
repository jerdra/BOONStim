nextflow.preview.dsl=2
import groovy.util.FileNameByRegexFinder


process project_mask2surf{

    label 'connectome'

    input:
    tuple val(sub), path(midthickness), path(white), path(pial), path(mask)

    output:
    tuple val(sub), path("surfmask.L.shape.gii"), emit: surfmask

    shell:
    '''
    #!/bin/bash

    wb_command -volume-to-surface-mapping \
                !{mask} \
                !{midthickness} \
                -ribbon-constrained \
                !{white} \
                !{pial} \
                "surfmask.L.shape.gii"
    '''

}

process clean_img{

    label 'ciftify'

    input:
    tuple val(sub), path(dtseries), path(confound), path(config)

    output:
    tuple val(sub), path("cleaned_img.dtseries.nii"), emit: clean_dtseries

    shell:
    '''
    ciftify_clean_img --clean-config=!{config} \
                        --confounds-tsv=!{confound} \
                        !{dtseries} \
                        --output-file cleaned_img.dtseries.nii
    '''
    
}

process smooth_img{

    label 'connectome'
        
    input:
    tuple val(sub), path(dtseries), path(left), path(right)
    
    output:
    tuple val(sub), path("smooth.dtseries.nii"), emit: smooth_dtseries

    shell:
    '''
    wb_command -cifti-smoothing \
                !{dtseries} \
                6 6 COLUMN \
                -left-surface !{left} \
                -right-surface !{right} \
                "smooth.dtseries.nii" 
    '''


}

process calculate_roi_correlation{

    label 'connectome'

    input:
    tuple val(sub), path(dtseries), path(shape)

    output:
    tuple val(sub), path("correlation.dscalar.nii"), emit: corr_dscalar

    shell:
    '''
    wb_command -cifti-average-roi-correlation \
                correlation.dscalar.nii \
                -left-roi !{shape} \
                -cifti !{dtseries}
    '''
}

process split_cifti{
    
    label 'connectome'
    
    input:
    tuple val(sub), path(dscalar)

    output:
    tuple val(sub), path('L.shape.gii'), emit: left_shape
    tuple val(sub), path('R.shape.gii'), emit: right_shape

    shell:
    '''
    wb_command -cifti-separate \
                !{dscalar} \
                COLUMN \
                -metric CORTEX_LEFT L.shape.gii \
                -metric CORTEX_RIGHT R.shape.gii
    '''


}

process mask_cortex{

    label 'connectome'
    
    input:
    tuple val(sub), path(shape)
    
    output:
    tuple val(sub), path("masked.${shape}"), emit: masked_shape

    shell:
    '''
    wb_command -metric-math \
                "x*0" \
                -var "x" !{shape} \
                masked.!{shape}
    '''

}

process create_dense{

    label 'connectome'
    
    input:
    tuple val(sub), path(left_shape), path(right_shape)

    output:
    tuple val(sub), path('weightfunc.dscalar.nii'), emit: weightfunc

    shell:
    '''
    wb_command -cifti-create-dense-scalar \
                weightfunc.dscalar.nii \
                -left-metric !{left_shape} \
                -right-metric !{right_shape}
    '''


}


workflow calculate_weightfunc_wf {

    get:
        derivatives

    main:
        
        // Project mask into surface space for the particular subject
        surfs = derivatives
                        .map{s,f,c ->   [
                                            s,
                                            "${s}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                            "${s}/MNINonLinear/fsaverage_LR32k/${s}.L.white.32k_fs_LR.surf.gii",
                                            "${s}/MNINonLinear/fsaverage_LR32k/${s}.L.pial.32k_fs_LR.surf.gii",
                                            "${params.inverse_mask}"
                                        ]
                            }
        project_mask2surf(surfs)

        // Get both dtseries files, split, get confounds and apply
        cleaned_input = derivatives
                            .map{s,f,c ->   [
                                                s,
                                                f,
                                                new FileNameByRegexFinder().getFileNames("${s}",".*MNINonLinear/Results/.*(REST|rest).*/.*dtseries.nii")
                                                
                                            ]
                                }
                            .transpose()
                            .map{ s,f,d ->  [
                                                s,f,d,
                                                ( d =~ /ses-[^_]*/ )[0],
                                                ( d =~ /ses-0._task.+?(?=(_desc|_Atlas))/ )[0][0]
                                            ]
                                }
                            .map{ s,f,d,ses,ident ->    [
                                                            s,d,
                                                            new FileNameByRegexFinder().getFileNames("$s/$ses/func", ".*${ident}.*confound.*tsv")[0],
                                                            "${params.clean_config}"
                                                        ]
                                }

        //Clean image
        clean_img(cleaned_input)
        
        //Smooth image need cifti information
        cifti_buffer = derivatives.map{s,f,c -> [s,c]}
        smooth_input = clean_img.out.clean_dtseries.join(cifti_buffer, by:0)
                                .map{ s,i,c ->  [
                                                    s,i,
                                                    "${s}/MNINonLinear/fsaverage_LR32k/${s}.L.midthickness.32k_fs_LR.surf.gii",
                                                    "${s}/MNINonLinear/fsaverage_LR32k/${s}.R.midthickness.32k_fs_LR.surf.gii"
                                                ]
                                    }
        smooth_img(smooth_input)

        //Compute correlation
        correlation_input = smooth_img.out.smooth_dtseries.join(project_mask2surf.out.surfmask, by:0 )
        calculate_roi_correlation(correlation_input)

        //Separate out cifti file
        split_cifti(calculate_roi_correlation.out.corr_dscalar)

        //Mask the right cortex entirely
        mask_cortex(split_cifti.out.right_shape)

        //Recombine to make final dscalar containing only left cortex
        create_dense_input = split_cifti.out.left_shape
                                        .join(mask_cortex.out.masked_shape, by: 0)
        create_dense(create_dense_input)
    
        emit:
            weightfunc = create_dense.out.weightfunc
        
}

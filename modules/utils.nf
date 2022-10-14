nextflow.preview.dsl=2

process apply_mask {

    /*
    Apply a mask to a dscalar image

    Arguments:
        sub (str): Subject ID
        dscalar (Path): Input dscalar file to mask
        mask (Path): Mask dscalar file

    Outputs:
        masked (channel): (sub, masked: Path)
    */

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(mask)

    output:
    tuple val(sub), path("${sub}.masked.dscalar.nii"), emit: masked

    shell:
    '''
    wb_command -cifti-math \
                "x * (mask > 0)" \
                -var "x" !{dscalar} \
                -var "mask" !{mask} \
                !{sub}.masked.dscalar.nii
    '''
}

process cifti_dilate {

    /*
    Dilate an input dscalar file

    Arguments:
        sub (str): Subject ID
        dscalar (Path): Path to dscalar file
        left (Path): Path to left surface file (midthickness typically)
        right (Path): Path to right surface file (midthickness typically)

    Outputs:
        dilated (channel): (sub, dilated: Path)
    */

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(left), path(right)

    output:
    tuple val(sub), path("${sub}.dilated.dscalar.nii"), emit: dilated

    shell:
    '''
    wb_command -cifti-dilate \
                !{dscalar} \
                COLUMN \
                6 6 \
                -left-surface !{left} \
                -right-surface !{right} \
                !{sub}.dilated.dscalar.nii
    '''
}

process numpy2txt {

    label 'fieldopt'

    input:
    tuple val(id), path(numpy)

    output:
    tuple val(id), path("numpyastxt.txt"), emit: txt

    shell:
    """
    #!/usr/bin/env python
    import numpy as np
    arr = np.load('${numpy}')
    np.savetxt('numpyastxt.txt', arr, delimiter=',')
    """
}

process txt2numpy {
    label 'fieldopt'

    input:
    tuple val(id), path(txt)

    output:
    tuple val(id), path("txtasnumpy.npy"), emit: npy

    shell:
    """
    #!/usr/bin/env python
    import numpy as np
    arr = np.load('${numpy}')
    arr.save('txtasnumpy.npy')
    """
}

process flip_direction_spec {
    label 'fieldopt'

    input:
    tuple val(id), path(json), path(flip)

    output:
    tuple val(id), path("${id}_flipped.json"), emit: spec

    shell:
    """
    #!/usr/bin/env python
    import json

    with open('${json}', 'r') as f:
        spec = json.load(f)

    with open('${flip}', 'r') as f:
        value = f.read()

    if value:
        value = int(value)

    if 'dir_x' in spec.keys() and value == 1:
        spec['dir_x'] = -spec['dir_x']
        spec['dir_y'] = -spec['dir_y']
        spec['dir_z'] = -spec['dir_z']

    with open('${id}_flipped.json', 'w') as f:
        json.dump(spec, f, indent=3)
    """
}

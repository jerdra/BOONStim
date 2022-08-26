import argparse

import numpy as np
from nilearn import image
import simnibs.msh.transformations as transformations

# SimNIBS convention for generating outputs
TMP_VOL_FILE = "fields.nii.gz"
TMP_NORME = "fields_normE.nii.gz"
TMP_TARGET = "fields_weightfunction.nii.gz"


def main():

    parser = argparse.ArgumentParser(
        description="Adjust dosage of E100 (100th largest E-field magnitude) "
        "to match that of a reference value")

    parser.add_argument("sim_msh", help="Simulation mesh", type=str)
    parser.add_argument("m2m_dir", help="mri2mesh directory", type=str)
    parser.add_argument("output", help="Output file", type=str)

    args = parser.parse_args()

    transformations.interpolate_to_volume(args.sim_msh, args.m2m_dir,
                                          TMP_VOL_FILE)

    roi = image.load_img(TMP_TARGET).get_fdata().astype(bool)
    roi_inds = np.where(roi)

    img = image.load_img(TMP_NORME)
    e_in_targ = img.get_fdata()[roi_inds]
    e100 = np.sort(e_in_targ)[-100]

    with open(args.output, 'w') as f:
        f.write(str(e100))


if __name__ == '__main__':
    main()

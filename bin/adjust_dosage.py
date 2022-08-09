import argparse
import json

import numpy as np
from nilearn import image
import simnibs.msh.transformations as transformations

# SimNIBS convention for generating outputs
TMP_VOL_FILE = "fields.nii.gz"
TMP_NORME = "fields_normE.nii.gz"
TMP_E = "fields_E.nii.gz"
TMP_TARGET = "fields_Target.nii.gz"


def main():

    parser = argparse.ArgumentParser(
        description="Adjust dosage of E100 (100th largest E-field magnitude) "
        "to match that of a reference value")

    parser.add_argument("sim_msh", help="Simulation mesh", type=str)
    parser.add_argument("m2m_dir", help="mri2mesh directory", type=str)
    parser.add_argument("reference", help="E100 reference value", type=float)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--direction-json",
                        help="Use a JSON file containing direction info",
                        type=str)

    args = parser.parse_args()

    transformations.interpolate_to_volume(args.sim_msh, args.m2m_dir,
                                          TMP_VOL_FILE)

    roi = image.load_img(TMP_TARGET).get_fdata().astype(bool)
    roi_inds = np.where(roi)

    if not args.direction_json:
        img = image.load_img(TMP_NORME)
        e_in_targ = img.get_fdata()[roi_inds]
        e100 = np.sort(e_in_targ)[-100]
    else:

        # Pull direction of radial
        with open(args.direction_json, 'r') as fhandle:
            spec = json.load(fhandle)
        direction = np.array([spec['dir_X'], spec['dir_Y'], spec['dir_Z']],
                             dtype=float).reshape((3, 1))

        # Ensure direction is unit vector
        direction = direction / np.linalg.norm(direction)

        # Get E-fields in target
        target_e = image.load_img(TMP_E).get_fdata()[roi_inds]

        # Compute magnitude of radial
        radial_mags = target_e @ direction
        e100 = np.sort(radial_mags)[-100]

    np.savetxt(args.output, e100)


if __name__ == '__main__':
    main()

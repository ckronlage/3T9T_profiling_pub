

import os
import glob
import numpy as np
import nibabel as nib

def main(outdir='data/derivatives/T2starmaps'):
    R2star_paths = glob.glob('data/sub-*/**/*R2starmap.nii.gz', recursive=True)

    for r_path in R2star_paths:
        target_path = f'{outdir}/{os.path.relpath(r_path, "data")}'.replace('R2starmap', 'T2starmap')
        print(f'Converting {r_path} to {target_path}')
        
        os.makedirs(os.path.dirname(target_path), exist_ok=True)

        r2star_img = nib.load(r_path)
        r2star_data = r2star_img.get_fdata()

        mask = r2star_data > 0
        t2star_data = np.zeros_like(r2star_data)
        t2star_data[mask] = 1000.0 / r2star_data[mask]  # T2* in ms where R2* > 0

        t2star_img = nib.Nifti1Image(t2star_data, affine=r2star_img.affine)
        nib.save(t2star_img, target_path)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert R2* maps to T2* maps.')
    parser.add_argument('--outdir', type=str, default='data/derivatives/T2starmaps',
                        help='Output directory for T2* maps')
    args = parser.parse_args()
    main(outdir=args.outdir)

import os
import argparse

import ants

import dataset_paths

def antsN4(filename_input,
           filename_output,
           iters=[50, 50, 50, 50]):
    print(f'N4 bias field correction, input: {filename_input}, iters: {iters}')
    if os.path.isfile(filename_output) and overwrite == False:
        print('already done, skipping')
        return
    image = ants.image_read(filename_input)
    # because ants N4 input cannot be negative and MP2RAGE is in the range [-0.5, 0.5], we rescale to [0, 1]:  
    image = image.iMath_normalize() + 0.01
    image_n4 = ants.n4_bias_field_correction(image, 
                                             convergence={'iters': iters, 'tol': 1e-5}, 
                                             verbose=False)
    ants.image_write(image_n4, filename_output)

    


def test_94T_N4_settings():
    outdir = 'data/derivatives/test_94T_N4_settings'
    # create output directories if they don't already exist
    os.makedirs(outdir, exist_ok=True)

    subjects = dataset_paths.get_3T94T_paths()
    subj = subjects[5]

    # first robust denoise input image
    BETA_94T = 0.4
    filename_robustMP2RAGE_94T = f'{outdir}/{subj.subj_id}_94T_robustMP2RAGE_{BETA_94T}.nii.gz'
    mp2rage_robust_combination(subj.UNIT1_94T,
                            subj.INV1_94T,
                            subj.INV2_94T,
                            filename_robustMP2RAGE_94T,
                            override_beta=BETA_94T,
                            overwrite=False)
    

    iters = [[50, 50, 50, 50],
             [200, 200, 200, 200],
             [500, 500, 500, 500],
             [50, 50, 50, 50, 50]]

    for i, params in enumerate(iters):
        filename_rob_MP2_N4_94T = f'{outdir}/{subj.subj_id}_94T_robustMP2RAGE_{BETA_94T}_N4_{i}.nii.gz'
        antsN4(filename_robustMP2RAGE_94T,
               filename_rob_MP2_N4_94T,
               iters=params)
        
    

            

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='ANTs N4 bias field correction (from antspyx)')
    parser.add_argument('filename_input', type=str, help='Input NIfTI file')
    parser.add_argument('filename_output', type=str, help='Output NIfTI file')
    parser.add_argument('--iters', type=int, nargs='+', default=[50, 50, 50, 50], help='N4 iterations')
    args = parser.parse_args()

    antsN4(args.filename_input, args.filename_output, iters=args.iters)

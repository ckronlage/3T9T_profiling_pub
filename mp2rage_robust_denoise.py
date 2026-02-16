
import argparse
from copy import deepcopy

import nibabel as nib
import numpy as np



def mp2rage_robust_combination(filename_uni, 
                               filename_inv1,
                               filename_inv2,
                               filename_output = None, 
                               multiplying_factor=1,
                               override_beta=None):
    """
    This is a version of the "robust" denoise MP2RAGE calculation described
    by O'Brien et al. (2014), originally implemented by Jose Marques in MATLAB
    ('RobustCombination.m', https://github.com/JosePMarques/MP2RAGE-related-scripts),
    ported to python by Marc Pabst (https://github.com/marcpabst/mp2rage-denoise/),
    now further adapted (keeping the original nifti header and estimating the noise-
    level at the most superior part of the image independent of axis orientation).

    Args:
        filename_uni (str): Path to the uniform T1-image (UNI).
        filename_inv1 (str): Path to the first inversion image (INV1).
        filename_inv2 (str): Path to the second inversion image (INV2).
        filename_output (str, optional): Path to output image.
        multiplying_factor (int, optional): Factor for calculating beta from the noise level 
                                            (as determined in an uppermost corner of the inv2 image).
        override_beta (float, optional): Fixed beta - if supplied, multiplying_factor is overriden
    """
    # define relevant functions
    mp2rage_robustfunc  =  lambda inv1, inv2, beta: (inv1.conj() * inv2 - beta) / (np.square(inv1) + np.square(inv2) + 2*beta)

    rootsquares_pos  = lambda a,b,c: (-b+np.sqrt(np.square(b) -4 *a*c))/(2*a)
    rootsquares_neg  = lambda a,b,c: (-b-np.sqrt(np.square(b) -4 *a*c))/(2*a)

    # load data
    image_uni  = nib.load(filename_uni)
    image_inv1 = nib.load(filename_inv1)
    image_inv2 = nib.load(filename_inv2)

    image_uni_fdata = image_uni.get_fdata()
    image_inv1_fdata = image_inv1.get_fdata()
    image_inv2_fdata  = image_inv2.get_fdata()

    # scale UNI image values 
    if (np.amin(image_uni_fdata) >=0) and (np.amax(image_uni_fdata >= 0.51)):
        scale = lambda x: (x - np.amax(image_uni_fdata)/2) / np.amax(image_uni_fdata)
        image_uni_fdata = scale(image_uni_fdata)
 
    # correct polarity for INV1
    image_inv1_fdata = np.sign(image_uni_fdata) * image_inv1_fdata

    # MP2RAGEimg is a phase sensitive coil combination.. some more maths has to
    # be performed to get a better INV1 estimate which here is done by assuming
    # both INV2 is closer to a real phase sensitive combination
    inv1_pos = rootsquares_pos(-image_uni_fdata, image_inv2_fdata, -np.square(image_inv2_fdata) * image_uni_fdata)
    inv1_neg = rootsquares_neg(-image_uni_fdata, image_inv2_fdata, -np.square(image_inv2_fdata) * image_uni_fdata)

    image_inv1_final_fdata = deepcopy(image_inv1_fdata)

    image_inv1_final_fdata[np.abs(image_inv1_fdata - inv1_pos) >  np.abs(image_inv1_fdata - inv1_neg)] = inv1_neg[np.abs(image_inv1_fdata - inv1_pos) >  np.abs(image_inv1_fdata - inv1_neg)]
    image_inv1_final_fdata[np.abs(image_inv1_fdata - inv1_pos) <= np.abs(image_inv1_fdata - inv1_neg)] = inv1_pos[np.abs(image_inv1_fdata - inv1_pos) <= np.abs(image_inv1_fdata - inv1_neg)]
    
    image_inv2_canonical = nib.as_closest_canonical(image_inv2) # so that the image is in RAS+ orientation
    noiselevel = multiplying_factor * np.mean(image_inv2_canonical.get_fdata()[:15,:20,-15:]) # one superior corner of the volume
    beta = np.square(noiselevel)
    if override_beta is not None:
        MAX12BIT = 4095.0
        beta = override_beta * MAX12BIT
    print(f'{filename_uni} - beta: {beta} - min: {image_inv2_fdata.min()} - max: {image_inv2_fdata.max()}')

    output = mp2rage_robustfunc(image_inv1_final_fdata, image_inv2_fdata, beta)
    image_output = nib.Nifti1Image(output, image_uni.affine, image_uni.header)

    if filename_output:
        nib.save(image_output, filename_output)
    else:
        return image_output
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MP2RAGE robust combination denoising')
    parser.add_argument('filename_uni', help='Path to the UNI image file')
    parser.add_argument('filename_inv1', help='Path to the INV1 image file')
    parser.add_argument('filename_inv2', help='Path to the INV2 image file')
    parser.add_argument('output_file', help='Path to the output file')
    parser.add_argument('--beta', type=float, default=0.4, help='beta value for the robust combination; if specified, overrides multiplying_factor')
    parser.add_argument('--multiplying_factor', type=float, default=1.0, help='parameter for estimating beta from the image noise level')
    args = parser.parse_args()

    mp2rage_robust_combination(args.filename_uni,
                               args.filename_inv1,
                               args.filename_inv2,
                               args.output_file,
                               multiplying_factor=args.multiplying_factor,
                               override_beta=args.beta)
    
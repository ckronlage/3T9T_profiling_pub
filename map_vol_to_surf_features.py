import os
import glob
import shutil
import argparse

import numpy as np
import scipy.ndimage 
import nibabel as nib
import ants

import dataset_paths
import util


def surface_features(freesurfer_dir, subj_id, threads=1):
    os.environ['SUBJECTS_DIR'] = freesurfer_dir

    vols_3T, vols_94T = coreg_other_vols(freesurfer_dir, subj_id, return_existing_vols=False, threads=threads)
    reestimate_thickness_max(freesurfer_dir, subj_id)
    reestimate_w_g_contrast(freesurfer_dir, f'sub-{subj_id}_3T')
    reestimate_w_g_contrast(freesurfer_dir, f'sub-{subj_id}_94T')
    sample_vol2surf_all_volumes(freesurfer_dir, subj_id, vols_3T, vols_94T)
    resample_surf_features_fsaverage_sym(freesurfer_dir,
                                            f'sub-{subj_id}_3T',
                                            f'{freesurfer_dir}/sub-{subj_id}_3T/new/surf')
    resample_surf_features_fsaverage_sym(freesurfer_dir,
                                            f'sub-{subj_id}_94T',
                                            f'{freesurfer_dir}/sub-{subj_id}_94T/new/surf')


def normalize_intensities(freesurfer_subj_folder, input_filename, output_filename):
    filename_aseg = f'{freesurfer_subj_folder}/mri/aseg.mgz'
    aseg = nib.load(filename_aseg)

    # select lh and rh white matter (indices 2 and 41) and cortex (indices 3 and 42)
    wm = np.isin(aseg.get_fdata(),[2, 41]) #.astype(float)
    cortex = np.isin(aseg.get_fdata(),[3, 42])

    # erode
    strel = scipy.ndimage.generate_binary_structure(3, 1)
    wm_eroded = scipy.ndimage.binary_erosion(wm, strel)
    cortex_eroded = scipy.ndimage.binary_erosion(cortex, strel)

    image = nib.load(input_filename)
    image_data = image.get_fdata()
    wm_median = np.median(image_data[wm_eroded])
    cortex_median = np.median(image_data[cortex_eroded])

    if 'FLAIR' in input_filename: # T2/FLAIR-contrast, i.e.: cortex has higher intensities than white matter
        normalized = (image_data - wm_median) / (cortex_median - wm_median)
    else:
        normalized = (image_data - cortex_median) / (wm_median - cortex_median)

    output = nib.Nifti1Image(normalized, image.affine, image.header)
    nib.save(output, output_filename)

def select_single_path(paths, subj_id, b0, sequence):
    path = paths[(paths['subj_id'] == subj_id) & (paths['B0'] == b0) & (paths['sequence'] == sequence)]['path'].values
    if len(path) == 1:
        return path[0]
    elif len(path) == 0:
        return None
    else:
        raise ValueError(f'More than 1 path found for {subj_id} {b0} {sequence}')

def coreg_other_vols(freesurfer_dir, subj_id, return_existing_vols=False, threads=1):
    paths = dataset_paths.get_3T94T_paths()

    # 3T #####################
    freesurfer_id_3T = f'sub-{subj_id}_3T'
    outdir = f'{freesurfer_dir}/{freesurfer_id_3T}/new'
    os.makedirs(outdir, exist_ok=True)

    vols_3T = []
    
    # FLAIR
    path_FLAIR_3T = select_single_path(paths, subj_id, '3T', 'FLAIR')
    outpath_FLAIR_3T = f'{outdir}/FLAIR_N4_to_orig_norm.nii.gz'
    if path_FLAIR_3T is not None:
        # N4 bias field correction
        FLAIR_N4_path = f'{outdir}/FLAIR_N4.nii.gz'
        im_FLAIR = ants.image_read(path_FLAIR_3T)
        im_FLAIR_N4 = ants.n4_bias_field_correction(im_FLAIR, 
                                            convergence={'iters': [50, 50, 50, 50], 'tol': 1e-7}, 
                                            verbose=False)
        ants.image_write(im_FLAIR_N4, FLAIR_N4_path) 
        
        # coreg with ref mask (because there is much extra-cerebral tissue in 
        # FLAIR and good grey-white contrast in the brain)
        FLAIR_reg_path = f'{outdir}/FLAIR_N4_to_orig_reg.lta'
        cmd = f'mri_coreg --mov {FLAIR_N4_path} --ref {outdir}/../mri/orig.mgz --reg {FLAIR_reg_path} --s {freesurfer_id_3T} --dof 12 --threads {threads}'
        util.bash_run(cmd)

        FLAIR_N4_to_orig_path = f'{outdir}/FLAIR_N4_to_orig.nii.gz'
        cmd = f'mri_vol2vol --mov {FLAIR_N4_path} --o {FLAIR_N4_to_orig_path} --reg {FLAIR_reg_path} --fstarg --trilinear'
        util.bash_run(cmd)

        # normalize
        normalize_intensities(f'{freesurfer_dir}/{freesurfer_id_3T}', FLAIR_N4_to_orig_path, outpath_FLAIR_3T)
        vols_3T.append(outpath_FLAIR_3T)

    # MPRAGE T1w
    path_T1w_3T = select_single_path(paths, subj_id, '3T', 'T1w')
    outpath_T1w_3T = f'{outdir}/T1w_N4_to_orig_norm.nii.gz'
    if path_T1w_3T is not None:
        # N4 bias field correction
        T1w_N4_path = f'{outdir}/T1w_N4.nii.gz'
        im_T1w = ants.image_read(path_T1w_3T)
        im_T1w_N4 = ants.n4_bias_field_correction(im_T1w, 
                                            convergence={'iters': [50, 50, 50, 50], 'tol': 1e-7}, 
                                            verbose=False)
        ants.image_write(im_T1w_N4, T1w_N4_path) 
        
        # coreg with ref mask (because there is much extra-cerebral tissue in
        # the T1w and good grey-white contrast in the brain)
        T1w_reg_path = f'{outdir}/T1w_N4_to_orig_reg.lta'
        cmd = f'mri_coreg --mov {T1w_N4_path} --ref {outdir}/../mri/orig.mgz --reg {T1w_reg_path} --s {freesurfer_id_3T} --dof 12 --threads {threads}'
        util.bash_run(cmd)

        T1w_N4_to_orig_path = f'{outdir}/T1w_N4_to_orig.nii.gz'
        cmd = f'mri_vol2vol --mov {T1w_N4_path} --o {T1w_N4_to_orig_path} --reg {T1w_reg_path} --fstarg --trilinear'
        util.bash_run(cmd)

        # normalize
        normalize_intensities(f'{freesurfer_dir}/{freesurfer_id_3T}', T1w_N4_to_orig_path, outpath_T1w_3T)
        vols_3T.append(outpath_T1w_3T)


    # MP2RAGE UNIT1 3T
    path_UNIT1_3T = select_single_path(paths, subj_id, '3T', 'UNIT1')
    outpath_UNIT1_3T = f'{outdir}/UNIT1_conform_orig.nii.gz'
    if path_UNIT1_3T is not None:
        # this needs to be output type float because mri_convert otherwise converts to uchar (8bit int)
        cmd = f'mri_convert {path_UNIT1_3T} {outdir}/UNIT1_conform_orig.nii.gz --conform -odt float'
        util.bash_run(cmd)
        vols_3T.append(outpath_UNIT1_3T)


    # T1map 3T
    path_T1map_3T = select_single_path(paths, subj_id, '3T', 'T1map')
    outpath_T1map_3T = f'{outdir}/T1map_conform_orig.nii.gz'
    if path_T1map_3T is not None:
        # this needs to be output type float because mri_convert otherwise converts to uchar (8bit int)
        cmd = f'mri_convert {path_T1map_3T} {outdir}/T1map_conform_orig.nii.gz --conform -odt float'
        util.bash_run(cmd)
        vols_3T.append(outpath_T1map_3T)

    
    # 94T #####################
    freesurfer_id_94T = f'sub-{subj_id}_94T'
    outdir = f'{freesurfer_dir}/{freesurfer_id_94T}/new'
    os.makedirs(outdir, exist_ok=True)

    vols_94T = []

    # MP2RAGE UNIT1 94T
    path_UNIT1_94T = select_single_path(paths, subj_id, '94T', 'UNIT1')
    outpath_UNIT1_94T = f'{outdir}/UNIT1_conform_orig.nii.gz'
    if path_UNIT1_94T is not None:
        # this needs to be output type float because mri_convert otherwise converts to uchar (8bit int)
        cmd = f'mri_convert {path_UNIT1_94T} {outdir}/UNIT1_conform_orig.nii.gz -cm -odt float'
        util.bash_run(cmd)
        vols_94T.append(outpath_UNIT1_94T)


    # T1map
    path_T1map_94T = select_single_path(paths, subj_id, '94T', 'T1map')
    outpath_T1map_94T = f'{outdir}/T1map_conform_orig.nii.gz'
    if path_T1map_94T is not None:
        # this needs to be output type float because mri_convert otherwise converts to uchar (8bit int)
        cmd = f'mri_convert {path_T1map_94T} {outdir}/T1map_conform_orig.nii.gz -cm -odt float'
        util.bash_run(cmd)
        vols_94T.append(outpath_T1map_94T)


    # QSM 3DFLASH echo 1, we need to coregister this to be able to align the QSM data
    path_GREecho1_94T = select_single_path(paths, subj_id, '94T', 'GREecho1')
    if path_GREecho1_94T is not None:
        # N4 bias field correction
        GREecho1_offline_N4_path = f'{outdir}/GREecho1_offline_N4.nii.gz'
        im_GREecho1_offline = ants.image_read(path_GREecho1_94T)
        im_GREecho1_offline_N4 = ants.n4_bias_field_correction(im_GREecho1_offline, 
                                            convergence={'iters': [50, 50, 50, 50], 'tol': 1e-7}, 
                                            verbose=False)
        ants.image_write(im_GREecho1_offline_N4, GREecho1_offline_N4_path) 
        
        # coreg without ref mask (because GRE is a slab only and not such a
        # good grey-white contrast in the brain), so skull and dura aid in 
        # registration
        GREecho1_offline_reg_path = f'{outdir}/GREecho1_offline_N4_to_orig_reg.lta'
        cmd = f'mri_coreg --mov {GREecho1_offline_N4_path} --ref {outdir}/../mri/orig.mgz --reg {GREecho1_offline_reg_path} --s {freesurfer_id_94T} --no-ref-mask --dof 12 --threads {threads}'
        util.bash_run(cmd)

        GREecho1_offline_N4_to_orig_path = f'{outdir}/GREecho1_offline_N4_to_orig.nii.gz'
        cmd = f'mri_vol2vol --mov {GREecho1_offline_N4_path} --o {GREecho1_offline_N4_to_orig_path} --reg {GREecho1_offline_reg_path} --fstarg --trilinear'
        util.bash_run(cmd)

        # QSM Tke 12
        path_QSMTke12_94T = select_single_path(paths, subj_id, '94T', 'QSMTke12')
        QSMTKe12_to_orig_path = f'{outdir}/QSMTke12_to_orig.nii.gz'
        cmd = f'mri_vol2vol --mov {path_QSMTke12_94T} --o {QSMTKe12_to_orig_path} --reg {GREecho1_offline_reg_path} --fstarg --trilinear'
        util.bash_run(cmd)
        vols_94T.append(QSMTKe12_to_orig_path)

        # QSM Tke 3
        path_QSMTke3_94T = select_single_path(paths, subj_id, '94T', 'QSMTke3')
        QSMTke3_to_orig_path = f'{outdir}/QSMTke3_to_orig.nii.gz'
        cmd = f'mri_vol2vol --mov {path_QSMTke3_94T} --o {QSMTke3_to_orig_path} --reg {GREecho1_offline_reg_path} --fstarg --trilinear'
        util.bash_run(cmd)
        vols_94T.append(QSMTke3_to_orig_path)

    # a different echo 1 in alignment with the R2star
    path_R2star_scanner_94T = select_single_path(paths, subj_id, '94T', 'R2star_scanner')
    path_GREecho1_scanner_94T = select_single_path(paths, subj_id, '94T', 'GREecho1_scanner')
    if path_R2star_scanner_94T is not None:
        if not subj_id == 'P009': # special case for P009, see below
            # coregister the scanner reconstruction echo 1 to the offline reconstruciton echo 1
            GREecho1_scanner_to_offline_reg_path = f'{outdir}/GREecho1_scanner_to_offline_reg.lta'
            cmd = f'mri_coreg --mov {path_GREecho1_scanner_94T} --ref {path_GREecho1_94T} --reg {GREecho1_scanner_to_offline_reg_path} --s {freesurfer_id_94T} --no-ref-mask --dof 12 --threads {threads}'
            util.bash_run(cmd)

            # combine the two registrations
            GREecho1_scanner_to_orig_reg_path = f'{outdir}/GREecho1_scanner_to_orig_reg.lta'
            cmd = f'mri_concatenate_lta {GREecho1_scanner_to_offline_reg_path} {GREecho1_offline_reg_path} {GREecho1_scanner_to_orig_reg_path}'
            util.bash_run(cmd)

            GREecho1_scanner_to_orig_path = f'{outdir}/GREecho1_scanner_to_orig.nii.gz'
            cmd = f'mri_vol2vol --mov {path_GREecho1_scanner_94T} --o {GREecho1_scanner_to_orig_path} --reg {GREecho1_scanner_to_orig_reg_path} --fstarg --trilinear'
            util.bash_run(cmd)

            # R2star is in register with the scanner reconstruction echo 1
            out = f'{outdir}/R2star_to_orig.nii.gz'
            cmd = f'mri_vol2vol --mov {path_R2star_scanner_94T} --o {out} --reg {GREecho1_scanner_to_orig_reg_path} --fstarg --trilinear'
            util.bash_run(cmd)
            vols_94T.append(out)

            # also the T2starmap
            path_T2star_scanner_94T = glob.glob(f'data/derivatives/T2starmaps/sub-{subj_id}/**/*T2starmap.nii.gz', recursive=True)
            if len(path_T2star_scanner_94T) == 1:
                path_T2star_scanner_94T = path_T2star_scanner_94T[0]
                out = f'{outdir}/T2star_to_orig.nii.gz'
                cmd = f'mri_vol2vol --mov {path_T2star_scanner_94T} --o {out} --reg {GREecho1_scanner_to_orig_reg_path} --fstarg --trilinear'
                util.bash_run(cmd)
                vols_94T.append(out)
            else:
                print(f'No unique T2star path found for subject {subj_id}')
        else: #special case for P009
        # for P009, the R2star was reconstructed offline, so is in register with
        # the GREecho1_94T and other files
            out = f'{outdir}/R2star_to_orig.nii.gz'
            cmd = f'mri_vol2vol --mov {path_R2star_scanner_94T} --o {out} --reg {GREecho1_offline_reg_path} --fstarg --trilinear'
            util.bash_run(cmd)
            vols_94T.append(out)

            # T2starmap
            path_T2star_scanner_94T = glob.glob(f'data/derivatives/T2starmaps/sub-{subj_id}/**/*T2starmap.nii.gz', recursive=True)
            if len(path_T2star_scanner_94T) == 1:
                path_T2star_scanner_94T = path_T2star_scanner_94T[0]
                out = f'{outdir}/T2star_to_orig.nii.gz'
                cmd = f'mri_vol2vol --mov {path_T2star_scanner_94T} --o {out} --reg {GREecho1_offline_reg_path} --fstarg --trilinear'
                util.bash_run(cmd)
                vols_94T.append(out)
            else:
                print(f'No unique T2star path found for subject {subj_id}')
    return vols_3T, vols_94T

def reestimate_thickness_max(freesurfer_dir, subj_id, max = 15.0):
    """resample cortical thickness with maximum of 15 mm (instead of default 5 mm)"""

    # 3T #####################
    freesurfer_id_3T = f'sub-{subj_id}_3T'
    outdir = f'{freesurfer_dir}/{freesurfer_id_3T}/new/surf'
    os.makedirs(outdir, exist_ok=True)
    
    for hemi in ['lh', 'rh']:
        if os.path.exists(f'{outdir}/nthickness_{hemi}.mgh'):
            continue
        cmd = f'mris_thickness -max {max} {freesurfer_id_3T} {hemi} {outdir}/nthickness_{hemi}.mgh'
        util.bash_run(cmd)
        # because mris_thickness always prepends a lh/rh to the output filename, rename the file
        os.rename(f'{outdir}/{hemi}.nthickness_{hemi}.mgh',f'{outdir}/nthickness_{hemi}.mgh')

    # 94T #####################
    freesurfer_id_94T = f'sub-{subj_id}_94T'
    outdir = f'{freesurfer_dir}/{freesurfer_id_94T}/new/surf'
    os.makedirs(outdir, exist_ok=True)

    for hemi in ['lh', 'rh']:
        if os.path.exists(f'{outdir}/nthickness_{hemi}.mgh'):
            continue
        cmd = f'mris_thickness -max {max} {freesurfer_id_94T} {hemi} {outdir}/nthickness_{hemi}.mgh'
        util.bash_run(cmd)
        # because mris_thickness always prepends a lh/rh to the output filename, rename the file
        os.rename(f'{outdir}/{hemi}.nthickness_{hemi}.mgh',f'{outdir}/nthickness_{hemi}.mgh')

def reestimate_w_g_contrast(freesurfer_dir, freesurfer_id):
    outdir = f'{freesurfer_dir}/{freesurfer_id}/new/surf'
    os.makedirs(outdir, exist_ok=True)
    
    UNIT1_nifti = f'{freesurfer_dir}/{freesurfer_id}/new/UNIT1_conform_orig.nii.gz'
    UNIT1_mgz = f'{freesurfer_dir}/{freesurfer_id}/new/UNIT1_conform_orig.mgz'
    if os.path.exists(UNIT1_nifti):
        UNIT1 = nib.load(UNIT1_nifti)
        UNIT1_data = UNIT1.get_fdata()
        if UNIT1_data.min() < 0: 
            # if the UNIT1 has negative values, it is probably in the range [-0.5,0.5]
            # needs to be rescaled to positive values
            UNIT1_data = UNIT1_data + 0.5
            # write back to new file
            UNIT1_nifti = f'{freesurfer_dir}/{freesurfer_id}/new/UNIT1_conform_orig_pos.nii.gz'
            nib.save(nib.Nifti1Image(UNIT1_data, UNIT1.affine, UNIT1.header), UNIT1_nifti)

        cmd = f'mri_convert {UNIT1_nifti} {UNIT1_mgz}'
        util.bash_run(cmd)

        cmd = f'pctsurfcon --s {freesurfer_id} --fsvol ../new/UNIT1_conform_orig --b w-g.pct.UNIT1_conform_orig'
        util.bash_run(cmd)

        for hemi in ['lh', 'rh']:
            surf_file = f'{freesurfer_dir}/{freesurfer_id}/surf/{hemi}.w-g.pct.UNIT1_conform_orig.mgh'
            shutil.move(surf_file,f'{outdir}/w-g.pct.UNIT1_{hemi}.mgh')
        


def sample_vol2surf(freesurfer_dir,
                    freesurfer_subjid,
                    source_vol,
                    projfracs,
                    projdists,
                    projabs):

    outdir = f'{freesurfer_dir}/{freesurfer_subjid}/new/surf'
    os.makedirs(outdir, exist_ok=True)
    source_base = os.path.basename(source_vol).replace('.nii.gz','')

    for hemi in ['lh', 'rh']:
        for projfrac in projfracs:
            out = f'{outdir}/{source_base}_projrel{projfrac}_{hemi}.mgh'
            if not os.path.isfile(out):
                cmd = f'mri_vol2surf --src {source_vol} --out {out} --cortex --hemi {hemi} --regheader {freesurfer_subjid} --projfrac {projfrac}'
                util.bash_run(cmd)

        for projdist in projdists:
            out = f'{outdir}/{source_base}_projrel{projdist}_{hemi}.mgh'
            if not os.path.isfile(out):
                cmd = f'mri_vol2surf --src {source_vol} --out {out} --cortex --hemi {hemi} --regheader {freesurfer_subjid} --projdist {projdist}'
                util.bash_run(cmd)

        for projab in projabs:
            out = f'{outdir}/{source_base}_projabs{projab}_{hemi}.mgh'
            if not os.path.isfile(out):
                cmd = f'mri_vol2surf --src {source_vol} --out {out} --cortex --hemi {hemi} --regheader {freesurfer_subjid} --projdist {projab} --surf pial'
                util.bash_run(cmd)

def sample_vol2surf_all_volumes(freesurfer_dir,
                                subj_id,
                                vols_3T,
                                vols_94T,
                                projfracs = [0.0, 0.25, 0.5, 0.75],
                                projdists = [-1.0, -2.0],
                                projabs = [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0]):
    for v in vols_3T:
        sample_vol2surf(freesurfer_dir,
                        f'sub-{subj_id}_3T',
                        v,
                        projfracs,
                        projdists,
                        projabs)
    
    for v in vols_94T:
        sample_vol2surf(freesurfer_dir,
                        f'sub-{subj_id}_94T',
                        v,
                        projfracs,
                        projdists,
                        projabs)


def resample_surf_features_fsaverage_sym(freesurfer_dir, freesufer_subjid, folder):
    surf_features_lh = glob.glob(f'{folder}/*lh.mgh')
    # filter out fsaverage_sym_lh files in case this is re-run
    surf_features_lh = [f for f in surf_features_lh if not 'fsaverage_sym' in f]

    surf_features_rh = glob.glob(f'{folder}/*rh.mgh')

    for f in surf_features_lh:
        base = f.replace('.mgh','')
        if os.path.isfile(f'{base}_fsaverage_sym_lh.mgh') and not 'abs' in base:
            continue
        cmd = f'mris_apply_reg --src {f} --trg {base}_fsaverage_sym_lh.mgh --streg {freesurfer_dir}/{freesufer_subjid}/surf/lh.fsaverage_sym.sphere.reg {freesurfer_dir}/fsaverage_sym/surf/lh.sphere.reg'
        util.bash_run(cmd)

    for f in surf_features_rh:
        base = f.replace('.mgh','')
        if os.path.isfile(f'{base}_fsaverage_sym_lh.mgh') and not 'abs' in base:
            continue
        cmd = f'mris_apply_reg --src {f} --trg {base}_fsaverage_sym_lh.mgh --streg {freesurfer_dir}/{freesufer_subjid}/xhemi/surf/lh.fsaverage_sym.sphere.reg {freesurfer_dir}/fsaverage_sym/surf/lh.sphere.reg'
        util.bash_run(cmd)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='map volume data to freesurfer surface reconstructions')
    argparser.add_argument('freesurfer_dir', type=str, help='Path to freesurfer directory')
    argparser.add_argument('subj_id', type=str, help='Subject ID (e.g., P001)')
    argparser.add_argument('--threads', type=int, default=1, help='Number of threads to use for coregistration (default: 1)')
    args = argparser.parse_args()

    surface_features(freesurfer_dir=args.freesurfer_dir,
                     subj_id=args.subj_id,
                     threads=args.threads)
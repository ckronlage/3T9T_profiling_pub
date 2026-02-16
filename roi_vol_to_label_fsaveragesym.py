
import os
import glob
import re
import argparse
import glob

import nibabel as nib
import numpy as np
import pygeodesic.geodesic as geodesic

import util


def roi_volume_to_fsaveragesym_label(freesurfer_subjid, roi_volume, outdir):
    hemis = ['lh', 'rh']
    for hemi in hemis:
        # map ROI from volume to surface
        roi_surf = f'{outdir}/{freesurfer_subjid}_{hemi}_roi_surf.mgh'
        cmd = f'mri_vol2surf --mov {roi_volume} --regheader {freesurfer_subjid} --hemi {hemi} --o {roi_surf} --projfrac-max -1 1 0.2'
        util.bash_run(cmd)

        # convert surface to label in native space
        roi_label = f'{outdir}/{freesurfer_subjid}_{hemi}_roi_label.mgh'
        cmd = f'mri_vol2label --i {roi_surf} --l {roi_label} --id 1 --surf {freesurfer_subjid} {hemi} white'
        util.bash_run(cmd)

    # map label to fsaverage_sym's left hemisphere
    path_label = glob.glob(f'{outdir}/{freesurfer_subjid}*roi_label.mgh.label')
    if not len(path_label) == 1:
        raise RuntimeError(f"Expecting one roi label, found: {path_label} for subject: {freesurfer_subjid}")
    path_label = path_label[0]

    hemi = re.sub(r'.*_(lh|rh)_.*',r'\1',path_label)

    if hemi == 'lh':
        path_label_fsaveragesym = f'{path_label}.lh_on_lh.fsaverage_sym.label'
        cmd = f'mris_apply_reg --src-label {path_label} --trg {path_label}.lh_on_lh.fsaverage_sym --streg $SUBJECTS_DIR/{freesurfer_subjid}/surf/lh.fsaverage_sym.sphere.reg $SUBJECTS_DIR/fsaverage_sym/surf/lh.sphere.reg '
    else:
        path_label_fsaveragesym = f'{path_label}.rh_on_lh.fsaverage_sym.label'
        cmd = f'mris_apply_reg --src-label {path_label} --trg {path_label}.rh_on_lh.fsaverage_sym --streg $SUBJECTS_DIR/{freesurfer_subjid}/xhemi/surf/lh.fsaverage_sym.sphere.reg $SUBJECTS_DIR/fsaverage_sym/surf/lh.sphere.reg '

    util.bash_run(cmd)
    return path_label_fsaveragesym, path_label


def align_roi_3T_to_94T(subj_id, roi_volume_3TMP2RAGE, denoised_3TMP2RAGE, outdir):
    freesurfer_subjid_3T = f'sub-{subj_id}_3T'
    freesurfer_subjid_94T = f'sub-{subj_id}_94T'

    # coregister denoised 3T MP2RAGE volume to 9.4T freesurfer orig volume
    reg = f'{outdir}/{freesurfer_subjid_3T}_to_94T.dat'
    if not os.path.isfile(reg): # because this takes longer
        cmd = f'bbregister --s {freesurfer_subjid_94T} --mov {denoised_3TMP2RAGE} --t1 --reg {reg}'
        util.bash_run(cmd)

    # apply registration to lesion ROI volume so that is aligned with sub-123_94T/mri/orig.mgz
    roi_volume_94TMP2RAGE = f'{outdir}/{freesurfer_subjid_94T}_roi_reg.mgz'
    cmd = f'mri_vol2vol --mov {roi_volume_3TMP2RAGE} --o {roi_volume_94TMP2RAGE} --reg {reg} --fstarg '
    util.bash_run(cmd)

    return roi_volume_94TMP2RAGE



def compute_lesion_distancemap_fsaveragesym(fs_base_subj_id, B0, hemi, path_lesionlabel):
    outdir = os.path.dirname(path_lesionlabel)
    path_distancemap = compute_lesion_distancemap(fs_base_subj_id, path_lesionlabel, hemi, B0)

    if hemi == 'lh':
        path_distancemap_fsaveragesym = f'{outdir}/{fs_base_subj_id}_{B0}_{hemi}_distances_lh_on_lh_fsaverage_sym.mgh'
        cmd = f'mris_apply_reg --src {path_distancemap} --trg {path_distancemap_fsaveragesym} --streg $SUBJECTS_DIR/{fs_base_subj_id}_{B0}/surf/lh.fsaverage_sym.sphere.reg $SUBJECTS_DIR/fsaverage_sym/surf/lh.sphere.reg '
        util.bash_run(cmd)
    else:
        path_distancemap_fsaveragesym = f'{outdir}/{fs_base_subj_id}_{B0}_{hemi}_distances_rh_on_lh_fsaverage_sym.mgh'
        cmd = f'mris_apply_reg --src {path_distancemap} --trg {path_distancemap_fsaveragesym} --streg $SUBJECTS_DIR/{fs_base_subj_id}_{B0}/xhemi/surf/lh.fsaverage_sym.sphere.reg $SUBJECTS_DIR/fsaverage_sym/surf/lh.sphere.reg '
        util.bash_run(cmd)

    return path_distancemap_fsaveragesym


def compute_lesion_distancemap(fs_base_subj_id, path_lesionlabel, hemi, B0):
    """
    Compute the geodesic distance map from the centroid of a lesion label on the cortical surface.

    Parameters:
    fs_base_subj_id (str): Subject identifier (sub-A123, without _B0 suffix)
    lesionlabel (str): Path to the lesion label file.
    hemi (str): Hemisphere ('lh' for left hemisphere, 'rh' for right hemisphere).
    B0 (str): Identifier for the B0 field strength.

    Returns:
    str: Path to the output file containing the distance map in MGH format.
    """

    # load surfaces
    pial_surf = nib.freesurfer.io.read_geometry(f'data/derivatives/freesurfer/{fs_base_subj_id}_{B0}/surf/{hemi}.pial')
    sphere_surf = nib.freesurfer.io.read_geometry(f'data/derivatives/freesurfer/{fs_base_subj_id}_{B0}/surf/{hemi}.sphere')

    # compute the distance from the centroid to the rest of the cortex
    lesionlabel = nib.freesurfer.io.read_label(path_lesionlabel)
    centroid_vertex = find_centroid(lesionlabel, sphere_surf)
    geoalg = geodesic.PyGeodesicAlgorithmExact(pial_surf[0], pial_surf[1])
    distances, _ = geoalg.geodesicDistances([centroid_vertex], list(range(len(pial_surf[0]))))

    # save to file
    outdir = os.path.dirname(path_lesionlabel)
    out_path = f'{outdir}/{fs_base_subj_id}_{B0}_{hemi}_distances.mgh'
    nib.save(nib.freesurfer.mghformat.MGHImage(distances.astype('float32'), affine=None), out_path)
    return out_path


def find_centroid(lesionlabel, sphere_surf):
    """
    Find the centroid of a lesion label on a spherical surface.

    Parameters:
    lesionlabel (array-like): Indices of the vertices that make up the lesion.
    sphere_ipsi (tuple): A tuple where the first element is an array of vertex coordinates.

    Returns:
    int: The index of the vertex that is the centroid of the lesion.
    """
    label_vertices = sphere_surf[0][lesionlabel]
    mean_coords = label_vertices.mean(axis=0)
    dists = [ (i, np.linalg.norm(sphere_surf[0][i] - mean_coords)) for i in lesionlabel ]
    dists_sorted = sorted(dists, key=lambda x: x[1])
    centroid = dists_sorted[0][0]
    return centroid


def roi_vol_to_label_fsaveragesym(subj_id, roi_volume_3TMP2RAGE, denoised_3TMP2RAGE, outdir = None):
    os.environ['SUBJECTS_DIR'] = os.path.abspath('data/derivatives/freesurfer')
    if outdir is None:
        outdir = f'{outdir}'
    os.makedirs(outdir,exist_ok=True)

    # map ROI volume in register with 3T MP2RAGE to 3T-based registration to fsaverage_sym
    path_label_3T_fsaveragesym, path_label_3T = roi_volume_to_fsaveragesym_label(freesurfer_subjid=f'sub-{subj_id}_3T',
                                                                                 roi_volume=roi_volume_3TMP2RAGE,
                                                                                 outdir=outdir)

    # align ROI volume to 94T MP2RAGE
    roi_vol_aligned_94T = align_roi_3T_to_94T(subj_id=subj_id,
                                                roi_volume_3TMP2RAGE=roi_volume_3TMP2RAGE,
                                                denoised_3TMP2RAGE=denoised_3TMP2RAGE,
                                                outdir=outdir)
    
    # map ROI volume in register with 94T MP2RAGE to 94T-based registration to fsaverage_sym
    path_label_94T_fsaveragesym, path_label_94T = roi_volume_to_fsaveragesym_label(freesurfer_subjid=f'sub-{subj_id}_94T',
                                                                                   roi_volume=roi_vol_aligned_94T,
                                                                                   outdir=outdir)

    # find hemi and output intersection between 3T and 94T versions
    hemi = re.sub(r'.*_(.+?)_roi_label.*',r'\1',path_label_94T_fsaveragesym)
    label_final_fsaverage_sym = f'{outdir}/sub-{subj_id}_{hemi}_final_fsaverage_sym_lh.label'
    cmd = f'labels_intersect {path_label_3T_fsaveragesym} {path_label_94T_fsaveragesym} {label_final_fsaverage_sym}'
    util.bash_run(cmd)

    #####################
    # compute distance maps

    path_distancemap_3T = compute_lesion_distancemap_fsaveragesym(f'sub-{subj_id}', '3T', hemi, path_label_3T)
    path_distancemap_94T = compute_lesion_distancemap_fsaveragesym(f'sub-{subj_id}', '94T', hemi, path_label_94T)

    # compute the mean of the two distance maps
    path_distancemap_mean = f'{outdir}/sub-{subj_id}_{hemi}_distance_mean.mgh'
    distancemap_mean = (nib.load(path_distancemap_3T).get_fdata() + nib.load(path_distancemap_94T).get_fdata()) / 2
    nib.save(nib.freesurfer.mghformat.MGHImage(distancemap_mean.astype('float32'), affine=None), path_distancemap_mean)


    #####################
    # write empty files for the other hemi
    empty_label = " #!ascii label \n 0 "
    other_hemi = 'lh' if hemi == 'rh' else 'rh'
    with open(f'{outdir}/sub-{subj_id}_{other_hemi}_final_fsaverage_sym_lh.label', 'w') as f:
        f.write(empty_label)

    path_distancemap_otherhemi = f'{outdir}/sub-{subj_id}_{other_hemi}_distance_mean.mgh'
    distancemap_nans = np.full_like(distancemap_mean, np.nan)
    nib.save(nib.freesurfer.mghformat.MGHImage(distancemap_nans.astype('float32'), affine=None), path_distancemap_otherhemi)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Convert ROI volume to label on fsaverage_sym and compute distance maps.')
    argparser.add_argument('subj_id', help='subject ID')
    argparser.add_argument('roi_volume_3TMP2RAGE', help='ROI volume in register with 3T MP2RAGE')
    argparser.add_argument('denoised_3TMP2RAGE', help='Denoised 3T MP2RAGE volume')
    argparser.add_argument('--outdir', default=None, help='Output directory')
    args = argparser.parse_args()

    roi_vol_to_label_fsaveragesym(subj_id=args.subj_id,
                                  roi_volume_3TMP2RAGE=args.roi_volume_3TMP2RAGE,
                                  denoised_3TMP2RAGE=args.denoised_3TMP2RAGE,
                                  outdir=args.outdir)

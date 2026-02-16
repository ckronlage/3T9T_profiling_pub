import glob
import os
import re


import nibabel as nib
import numpy as np
import pandas as pd



def read_surfdata_and_lesions():
    lesions = glob.glob(f'data/derivatives/lesionlabels/sub*/*final_fsaverage_sym_lh.label')
    lesional_subjects = [re.search(r'sub-P\d{3}',lesion)[0] for lesion in lesions]
    lesional_subjects = list(set(lesional_subjects)) # unique subject_ids

    lesional_freesurfer_dirs = [glob.glob(f'data/derivatives/freesurfer/{subject}*') for subject in lesional_subjects]
    lesional_freesurfer_dirs = sum(lesional_freesurfer_dirs, []) # flatten the list
    control_freesurfer_dirs = glob.glob(f'data/derivatives/freesurfer/sub-C*')

    surfdata = read_all_surfdata(lesional_freesurfer_dirs + control_freesurfer_dirs)

    # read lesion labels
    for lesion in lesions:
        indices = nib.freesurfer.io.read_label(lesion)
        if len(indices) == 0:
            continue
        boolean = np.zeros(surfdata.shape[1],dtype=bool)
        boolean[indices] = True

        subject = re.search(r'sub-P\d{3}',lesion)[0]
        hemi = re.search(r'(lh|rh)',lesion)[0]

        surfdata.loc[(subject, np.nan, 'lesion_label', hemi, np.nan)] = boolean

    # read distance maps
    distance_maps = glob.glob(f'data/derivatives/lesionlabels/sub-*/*distance_mean.mgh')
    for distance_map in distance_maps:
        distance_data = nib.load(distance_map).get_fdata()[:,0,0]
        if np.all(np.isnan(distance_data)):
            continue
        subject = re.search(r'sub-P\d{3}',distance_map)[0]
        hemi = re.search(r'(lh|rh)',distance_map)[0]

        surfdata.loc[(subject, np.nan, 'lesion_pial_distance', hemi, np.nan)] = distance_data

    # read fsaverage_sym cortex label
    indices = nib.freesurfer.io.read_label('data/derivatives/freesurfer/fsaverage_sym/label/lh.cortex.label')
    boolean = np.zeros(surfdata.shape[1],dtype=bool)
    boolean[indices] = True
    surfdata.loc[('fsaverage_sym', np.nan, 'cortex_label', np.nan, np.nan)] = boolean

    return surfdata

def read_all_surfdata(list_of_freesurfer_dirs):
    metadata_list = []
    data_list = []
    for control in list_of_freesurfer_dirs:
        print(control)
        data, metadata = read_surf_features(path_to_fs_subj_dir=control)
        data_list.extend(data)
        metadata_list.extend(metadata)
    data_ndarray = np.array(data_list)
    index = pd.MultiIndex.from_tuples(metadata_list, names=['subject', 'B0', 'feature', 'hemi', 'depth'])
    dataframe = pd.DataFrame(data_ndarray, index=index)

    # for QSM, R2star and T2star in controls, replace data are 0.0 (not contained in the
    # measured slab) with NaN
    i = dataframe.index.to_frame()
    sel = (i['subject'].str.contains('C')) & (i['feature'].str.contains('QSMTke3|QSMTke12|R2star|T2star'))
    dataframe.loc[sel] = dataframe.loc[sel].replace(0.0, np.nan)

    # rescale 3T T1map from ms to s
    i = dataframe.index.to_frame()
    sel = (i['B0'] == '3T') & (i['feature'].str.contains('T1map'))
    dataframe.loc[sel] = dataframe.loc[sel] / 1000.0

    return dataframe

def read_surf_features(path_to_fs_subj_dir):
    metadata_list = []
    data_list = []
    subject = re.search(r'sub-(C|P)\d{3}',os.path.basename(path_to_fs_subj_dir))[0]
    B0 = re.search(r'(3T|94T)',os.path.basename(path_to_fs_subj_dir))[0]

    features = glob.glob(f'data/derivatives/freesurfer/{subject}_{B0}/new/surf/*_fsaverage_sym_lh.mgh')

    for feature in features:
        feature_name = re.sub(r'(.*?)_.*(?:lh|rh).*',r'\1',os.path.basename(feature))
        hemi = re.search(r'(lh|rh)',os.path.basename(feature))[0]

        if 'proj' in feature:
            feature_depth = re.sub(r'.*proj(?:rel|abs)(.*?)_.*',r'\1',os.path.basename(feature))
            feature_depth = float(feature_depth)
            if 'projabs' in feature:
                feature_name += '_projabs'
            if 'projrel' in feature:
                feature_name += '_projrel'
        else:
            feature_depth = np.nan

        data = nib.load(feature).get_fdata()
        if data.shape[1] != 1 or data.shape[2] != 1:
            raise RuntimeError(f"Expected surface feature file, found wrong dimensions of {data.shape} in file {feature}")
                
        surf_data = data[:,0,0]
        data_list.append(surf_data)
        
        metadata = [subject, B0, feature_name, hemi, feature_depth]
        metadata_list.append(metadata)

    return data_list, metadata_list


if __name__ == '__main__':
    surfdata = read_surfdata_and_lesions()
    surfdata.to_pickle('tmp/surfdata.pkl')
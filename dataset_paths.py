import os
import glob
import re
import pandas as pd



def get_3T94T_paths(datapath='data'):
    """
    Retrieve paths for 3T and 9.4T MRI data for subjects in a given dataset directory.

    """

    if not os.path.isdir(datapath):
        raise ValueError(f'{datapath} is not an existing folder')

    # filter unique subject IDs from filenames using regex
    subj_ids = []
    folders = glob.glob(os.path.join(datapath, 'sub*'))
    for f in folders:
        subj_id = re.search(r'.*sub-(....).*', os.path.basename(f))
        if subj_id is not None and subj_id.group(1) not in subj_ids:
            subj_ids.append(subj_id.group(1))

    # check whether each subject has a ses-3T and ses-94T subfolder
    subj_ids = [s for s in subj_ids 
                if (os.path.isdir(os.path.join(datapath, f'sub-{s}', 'ses-3T')) 
                    and 
                    os.path.isdir(os.path.join(datapath, f'sub-{s}', 'ses-94T')))]

    all_nifti_paths = glob.glob(f'{datapath}/sub-*/**/*.nii*', recursive=True)

    paths = []
    # look for relevant paths for each subject
    for subj_id in subj_ids:
        # filter out relevant paths first (speeds things up a lot)
        regex = re.compile(f'sub-{subj_id}')
        subj_nifti_paths = [p for p in all_nifti_paths if re.search(regex, p)]

        def filter_paths(subj_id, B0, sequence, pattern):
            regex = re.compile(pattern)
            matches = [p for p in subj_nifti_paths if re.search(regex, p)]
            if len(matches) > 1:
                raise ValueError(f'More than 1 file found for pattern: {pattern}')
            if len(matches) == 0:
                return None
            return [subj_id, B0, sequence, matches[0]]

        paths.append(filter_paths(subj_id, '3T', 'UNIT1', f'sub-{subj_id}.*ses-3T.*anat.*acq-A_UNIT1.*.nii'))
        paths.append(filter_paths(subj_id, '3T', 'INV1', f'sub-{subj_id}.*ses-3T.*anat.*acq-A_inv-1.*.nii'))
        paths.append(filter_paths(subj_id, '3T', 'INV2', f'sub-{subj_id}.*ses-3T.*anat.*acq-A_inv-2.*.nii'))
        paths.append(filter_paths(subj_id, '3T', 'T1map', f'sub-{subj_id}.*ses-3T.*anat.*acq-A_T1map.nii'))
        paths.append(filter_paths(subj_id, '3T', 'T1w', f'sub-{subj_id}.*ses-3T.*anat.*acq-MPRAGE_T1w.*.nii'))
        paths.append(filter_paths(subj_id, '3T', 'FLAIR', f'sub-{subj_id}.*ses-3T.*anat.*acq-T2SPACE_FLAIR.*.nii'))
        paths.append(filter_paths(subj_id, '3T', 'ROI', f'sub-{subj_id}.*ses-3T.*anat.*roi.*.nii'))

        paths.append(filter_paths(subj_id, '94T', 'UNIT1', f'sub-{subj_id}.*ses-94T.*anat.*acq-old_UNIT1.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'INV1', f'sub-{subj_id}.*ses-94T.*anat.*acq-old_inv-1.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'INV2', f'sub-{subj_id}.*ses-94T.*anat.*acq-old_inv-2.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'T1map', f'sub-{subj_id}.*ses-94T.*anat.*acq-old_T1map.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'QSMTke3', f'sub-{subj_id}.*ses-94T.*anat.*acq-3DFLASH_rec-Tke3_proc-offline.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'QSMTke12', f'sub-{subj_id}.*ses-94T.*anat.*acq-3DFLASH_rec-Tke12_proc-offline.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'GREecho1', f'sub-{subj_id}.*ses-94T.*anat.*acq-3DFLASH_echo-1_proc-offline_T2starw.*.nii'))
        paths.append(filter_paths(subj_id, '94T', 'R2star_scanner', f'sub-{subj_id}.*ses-94T.*anat.*acq-3DFLASH.*R2starmap.*'))
        paths.append(filter_paths(subj_id, '94T', 'GREecho1_scanner', f'sub-{subj_id}.*ses-94T.*anat.*acq-3DFLASH_echo-1_proc-scanner_T2starw.*.nii'))

    paths = [p for p in paths if p is not None]
    paths = pd.DataFrame(paths, columns=['subj_id', 'B0', 'sequence', 'path'])

    return paths

def main():
    paths = get_3T94T_paths()
    pass

if __name__ == '__main__':
    main()

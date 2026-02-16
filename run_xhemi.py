import os
import glob
import util


def freesurfer_xhemi(dry_run : bool,
                    qsub : bool,
                    freesurfer_subjects_dir : str, 
                    freesurfer_subj_id : str):
    
    freesurfer_subjects_dir = os.path.abspath(freesurfer_subjects_dir)

    if not os.path.isdir(f'{freesurfer_subjects_dir}/{freesurfer_subj_id}'):
        raise ValueError(f"freesurfer subject folder not found, looking for: {freesurfer_subjects_dir}/{freesurfer_subj_id}")

    if os.path.isdir(f'{freesurfer_subjects_dir}/{freesurfer_subj_id}/xhemi'):
        print(f"xhemi apparently already run for subject {freesurfer_subj_id} - skipping")
        return

    script_cmd = f'''
        #!/bin/bash
        source nic_env.sh
        SetFreeSurfer 7.3.2
        export SUBJECTS_DIR={freesurfer_subjects_dir}
        cd {freesurfer_subjects_dir}
        surfreg --s {freesurfer_subj_id} --t fsaverage_sym --lh --no-annot
        surfreg --s {freesurfer_subj_id} --t fsaverage_sym --lh --no-annot --xhemi
        '''
    
    logdir = f'{freesurfer_subjects_dir}/logs/{freesurfer_subj_id}/'
    os.makedirs(logdir, exist_ok=True)

    util.write_script_to_file(f'{logdir}/xhemi.sh', script_cmd)

    cmd = f'{logdir}/xhemi.sh'
    if qsub:    
        cmd = f'''qsub 
                -q long.q 
                -cwd -V -b y
                -e {logdir}
                -o {logdir}
                {cmd}'''
        cmd = cmd.replace('\n',' ')

    if dry_run:
        print(cmd)
    else:
        util.bash_run(cmd)


def run_all_freesurfer_xhemi(freesurfer_subjects_dir, dry_run=False, qsub=True):
    """
    Run xhemi surface reconstructions, this needs a  separate script because of a little
    typo in surfreg in 7.4.1 (https://github.com/freesurfer/freesurfer/commit/69eb32330500e4f566f6ff0bae119dc3e9ead037),
    so we use 7.3.2

    """

    fs_subjs = glob.glob(f'{freesurfer_subjects_dir}/sub-*')
    fs_subjs = [os.path.split(p)[1] for p in fs_subjs] 

    for fs_sub in fs_subjs:
        freesurfer_xhemi(dry_run=dry_run,
                         qsub=qsub,
                         freesurfer_subjects_dir=freesurfer_subjects_dir,
                         freesurfer_subj_id=fs_sub)
        
def main():
    freesurfer_subjects_dir = os.path.abspath('data/derivatives/freesurfer')
    run_all_freesurfer_xhemi(freesurfer_subjects_dir)

if __name__=='__main__':
    main()

import os


import util


def main():
    #subj_id = 'P026'
    #center = (-21.43, 52.40, 23.38)
    #plot_lesion_3T94T(subj_id, center)

    subj_id = 'P017'
    center = (18.23, 54.17, 20.78)
    plot_lesion_3T94T(subj_id, center)

    

def plot_lesion_3T94T(subj_id, center):
    # Create tmp directory if it doesn't exist
    os.makedirs('tmp', exist_ok=True)

    os.environ['SUBJECTS_DIR'] = 'data/derivatives/freesurfer'

    dir_3T = f'data/derivatives/freesurfer/sub-{subj_id}_3T/'
    dir_94T = f'data/derivatives/freesurfer/sub-{subj_id}_94T/'

    path_coreg_3T_to_94T = f'data/derivatives/freesurfer/sub-{subj_id}_3T/new/3T_to_94T.lta'
    if not os.path.exists(path_coreg_3T_to_94T):
        cmd = f'bbregister --s sub-{subj_id}_94T --mov data/derivatives/freesurfer/sub-{subj_id}_3T/mri/orig.mgz --reg data/derivatives/freesurfer/sub-{subj_id}_3T/new/3T_to_94T.lta --t1 --12'
        util.bash_run(cmd)

    cmd = rf'''
        freeview -v \
            {dir_94T}/new/UNIT1_conform_orig.mgz:name=94T_UNIT1 \
            {dir_3T}/new/UNIT1_conform_orig.mgz:reg={path_coreg_3T_to_94T}:resample=cubic:name=3T_UNIT1 \
            --ras {center[0]} {center[1]} {center[2]} \
        -f \
            {dir_94T}/surf/lh.white:edgecolor=blue \
            {dir_94T}/surf/rh.white:edgecolor=blue \
            {dir_94T}/surf/lh.pial:edgecolor=cyan \
            {dir_94T}/surf/rh.pial:edgecolor=cyan \
            {dir_3T}/surf/lh.white:affinexfm={path_coreg_3T_to_94T}:edgecolor=blue \
            {dir_3T}/surf/rh.white:affinexfm={path_coreg_3T_to_94T}:edgecolor=blue \
            {dir_3T}/surf/lh.pial:affinexfm={path_coreg_3T_to_94T}:edgecolor=cyan \
            {dir_3T}/surf/rh.pial:affinexfm={path_coreg_3T_to_94T}:edgecolor=cyan \
            --layout 4
    '''
    print(cmd)
    util.bash_run(cmd)


if __name__ == '__main__':
    main()
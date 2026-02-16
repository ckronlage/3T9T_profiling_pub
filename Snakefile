import os
 
# allows changing the shell_prefix from the workflow profile
if 'shell_prefix' in config.keys():
    shell.prefix(config['shell_prefix'])

import dataset_paths
subjects = dataset_paths.get_3T94T_paths()


wildcard_constraints:
    # no / or _ in subj_id
    subj_id=r"[^/_]+"

rule mp2rage_robust_denoise:
    input:
        uni = ancient(lambda wildcards: subjects[(subjects['subj_id'] == wildcards.subj_id) & 
                                   (subjects['B0'] == wildcards.B0) &
                                   (subjects['sequence'] == 'UNIT1')]
                                   ['path'].squeeze()),
        inv1 = ancient(lambda wildcards: subjects[(subjects['subj_id'] == wildcards.subj_id) & 
                                   (subjects['B0'] == wildcards.B0) &
                                   (subjects['sequence'] == 'INV1')]
                                   ['path'].squeeze()),
        inv2 = ancient(lambda wildcards: subjects[(subjects['subj_id'] == wildcards.subj_id) & 
                                   (subjects['B0'] == wildcards.B0) &
                                   (subjects['sequence'] == 'INV2')]
                                   ['path'].squeeze()),
    output:
        '{path}/{subj_id}_{B0}_robustMP2RAGE_{beta}.nii.gz'
    wildcard_constraints:
        # no underscore in beta
        beta=r"[^/_]+"
    conda:
        '3T9T_profiling'
    localrule:
        True
    shell:
        """
        python mp2rage_robust_denoise.py {input.uni} {input.inv1} {input.inv2} {output} --beta {wildcards.beta}
        """

rule antsN4:
    input:
        'data/derivatives/proc/{subj_id}_{B0}_robustMP2RAGE_0.4.nii.gz'
    output:
        'data/derivatives/proc/{subj_id}_{B0}_robustMP2RAGE_0.4_N4.nii.gz'
    conda:
        '3T9T_profiling'
    params:
        iters='200 200 200 200'
    shell:
        """
        python antsN4.py {input} {output} --iters {params.iters}
        """

rule recon_all_3T:
    input:
        'data/derivatives/proc/{subj_id}_3T_robustMP2RAGE_0.4_N4.nii.gz'
    output:
        done='data/derivatives/freesurfer/sub-{subj_id}_3T/scripts/recon-all.done',
        aseg='data/derivatives/freesurfer/sub-{subj_id}_3T/mri/aseg.mgz'
    threads:
        4
    container:
        "docker://freesurfer/freesurfer:7.4.1"
    shell:
        """
        rm -rf data/derivatives/freesurfer/sub-{wildcards.subj_id}_3T
        export FS_ALLOW_DEEP=1
        export SUBJECTS_DIR=data/derivatives/freesurfer
        echo "synthstrip --no-csf" > data/derivatives/freesurfer/expert_3T.opts
        recon-all -s sub-{wildcards.subj_id}_3T -sd data/derivatives/freesurfer -i {input} -all -synthstrip -expert data/derivatives/freesurfer/expert_3T.opts -openmp {threads}
        """

rule recon_all_94T:
    input:
        'data/derivatives/proc/{subj_id}_94T_robustMP2RAGE_0.4_N4.nii.gz'
    output:
        done='data/derivatives/freesurfer/sub-{subj_id}_94T/scripts/recon-all.done',
        aseg='data/derivatives/freesurfer/sub-{subj_id}_94T/mri/aseg.mgz'
    threads:
        4
    container:
        "docker://freesurfer/freesurfer:7.4.1"
    shell:
        """
        rm -rf data/derivatives/freesurfer/sub-{wildcards.subj_id}_94T
        export FS_ALLOW_DEEP=1
        export SUBJECTS_DIR=data/derivatives/freesurfer
        echo "synthstrip --no-csf" > data/derivatives/freesurfer/expert_94T.opts
        echo "mris_inflate -n 50" >> data/derivatives/freesurfer/expert_94T.opts
        recon-all -s sub-{wildcards.subj_id}_94T -sd data/derivatives/freesurfer -i {input} -all -synthstrip -expert data/derivatives/freesurfer/expert_94T.opts -openmp {threads} -hires
        """

rule xhemi:
    input:
        'data/derivatives/freesurfer/sub-{subj_id}_{B0}/scripts/recon-all.done'
    output:
        'data/derivatives/freesurfer/sub-{subj_id}_{B0}/surf/lh.fsaverage_sym.sphere.reg',
        'data/derivatives/freesurfer/sub-{subj_id}_{B0}/xhemi/surf/lh.fsaverage_sym.sphere.reg'
    container:
        "docker://freesurfer/freesurfer:7.3.2"
    shell:
        """
        export SUBJECTS_DIR=$(realpath data/derivatives/freesurfer)
        cd $SUBJECTS_DIR
        rm -r sub-{wildcards.subj_id}_{wildcards.B0}/xhemi
        surfreg --s sub-{wildcards.subj_id}_{wildcards.B0} --t fsaverage_sym --lh --no-annot
        surfreg --s sub-{wildcards.subj_id}_{wildcards.B0} --t fsaverage_sym --lh --no-annot --xhemi
        """

rule lesion_roi_vol_to_label_fsaveragesym:
    input:
        roi_volume_3TMP2RAGE = ancient(lambda wildcards: subjects[(subjects['subj_id'] == wildcards.subj_id) & 
                                                                  (subjects['B0'] == '3T') &
                                                                  (subjects['sequence'] == 'ROI')]
                                                                  ['path'].squeeze()),
        denoised_3TMP2RAGE = 'data/derivatives/proc/{subj_id}_3T_robustMP2RAGE_0.4_N4.nii.gz'
    output:
        'data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_lh_final_fsaverage_sym_lh.label',
        'data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_rh_final_fsaverage_sym_lh.label',
        'data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_lh_distance_mean.mgh',
        'data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_rh_distance_mean.mgh'
    conda:
        '3T9T_profiling'
    shell:
        """
        rm -r data/derivatives/lesionlabels/sub-{wildcards.subj_id}/
        python roi_vol_to_label_fsaveragesym.py {wildcards.subj_id} {input.roi_volume_3TMP2RAGE} {input.denoised_3TMP2RAGE} --outdir data/derivatives/lesionlabels/sub-{wildcards.subj_id}/
        """

rule r2star_to_t2starmaps:
    input:

    output:
        directory('data/derivatives/T2starmaps/')
    conda:
        '3T9T_profiling'
    shell:
        """
        python R2star_to_T2star.py --outdir {output}
        """

rule map_vol_to_surf_features:
    input:
        'data/derivatives/freesurfer/sub-{subj_id}_{B0}/scripts/recon-all.done',
        'data/derivatives/T2starmaps/'
    output:
        'data/derivatives/freesurfer/sub-{subj_id}_{B0}/new/surf/UNIT1_conform_orig_projrel0.0_lh_fsaverage_sym_lh.mgh' # just one of the outputs that exists for all subjects
    conda:
        '3T9T_profiling'
    shell:
        """
        rm -r data/derivatives/freesurfer/sub-{wildcards.subj_id}_{wildcards.B0}/new/
        python map_vol_to_surf_features.py data/derivatives/freesurfer {wildcards.subj_id}
        """

rule gather_surf_features:
    input:
        # features vol to surf
        expand('data/derivatives/freesurfer/sub-{subj_id}_{B0}/new/surf/UNIT1_conform_orig_projrel0.0_lh_fsaverage_sym_lh.mgh',
                subj_id = subjects['subj_id'].unique(),
                B0 = ['3T', '94T']),

        # lesion ROIs
        expand('data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_{hemi}_final_fsaverage_sym_lh.label',
                subj_id = subjects[subjects['sequence']=='ROI']['subj_id'].unique(),
                hemi = ['lh', 'rh']),
        expand('data/derivatives/lesionlabels/sub-{subj_id}/sub-{subj_id}_{hemi}_distance_mean.mgh',
                subj_id = subjects[subjects['sequence']=='ROI']['subj_id'].unique(),
                hemi = ['lh', 'rh'])      
    output:
        'tmp/surfdata.pkl'
    conda:
        '3T9T_profiling'
    shell:
        """
        python gather_surf_features.py
        """

rule all:
    input:
        'tmp/surfdata.pkl'
    message:
        "rule {rule}"
    localrule:
        True
    default_target:
        True


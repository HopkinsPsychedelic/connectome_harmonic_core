# CHAP (Connectome Harmonic Analysis Pipeline)
CHAP (Connectome Harmonic Analysis Pipeline) is a containerized software pipeline for performing harmonic decomposition of structural connectomes. CHAP contains two separate arms: CHAP-HCP and CHAP-BIDS, which refer to the type of data that each arm accepts. 

CHAP-HCP takes as input Human Connectome Project data preprocessed according to the HCP minimal preprocessing pipeline. For each set of data from each subject and session, data contained in the HCP Structural, Diffusion, and Structural Extended downloads are used to compute structural harmonics. HCP resting state and task-based fMRI data are optional; when present, CHAP-HCP will output a variety of functional spectra for these scans. 

CHAP-BIDS can compute connectome harmonics on any dataset containing T1w and diffusion data with reverse phase-encoding fieldmaps. At this time, this arm requires the outputs of three open-source BIDS (Brain Imaging Data Structure) apps: *mrtrix3_connectome*, for preprocessing diffusion data,  *fmriprep* for preprocessing of fMRI data and running FreeSurfer, and *ciftify*, for further processing of surfaces and functional data. 

**Notes on fMRIprep**:

[fMRIprep](https://fmriprep.org/en/stable/) should be run first. Make sure the "anat" output space is selected, and make sure FreeSurfer is run.

**Notes on ciftify**:

[Ciftify](https://edickie.github.io/ciftify/#/) should be run with the following two options: 1) --read-from-derivatives should point to the derivatives directory that contains your fMRIprep and FreeSurfer output folders. 2)--resample-to-T1w32k 

**Notes on MRtrix3_connectome**:

[MRtrix3_connectome](https://github.com/BIDS-Apps/MRtrix3_connectome) should be run with the following options: 1) the "preproc" level should be run as opposed to the "participant" level. 2) --t1w_preproc should point to the location of the desc-preproc T1w in the fMRIprep output directory. 

**Running CHAP**

CHAP can be run through Docker or Singularity. To use singularity, create a .sif using: singularity pull docker://winstonian3/connectome_harmonic

For good information on running Docker and Singularity apps, [see this](https://www.nipreps.org/apps/docker/#running-a-niprep-directly-interacting-with-the-docker-engine).

In the below example, I run CHAP-HCP. The two required positional arguments are output directory and analysis level, respectively. For now, only "participant" level is supported. In CHAP-HCP, the --hcp_dir argument is also required. If using test-retest data, the hcp_dir should contain two directories: ses-test and ses-retest, each with the data directories downloaded as-is from ConnectomeDB. If just one session, the HCP downloads can go directly in hcp_dir. The --participant_label argument is optional. If this parameter is not provided all subjects will be analyzed. 

    sudo docker run -it --rm \
    -v /data/hcp_test_retest_pp/:/data/hcp_test_retest_pp/  \
    -p 8888:8888 \
    winstonian3/connectome_harmonic:latest \
    /data/hcp_test_retest_pp/derivatives \
    participant \
    --hcp_dir /data/hcp_test_retest_pp/source_data \
    --participant_label 105923

In the below example, I run CHAP-BIDS. In addition to the two positional arguments discussed above, --mrtrix_dir, --ciftify_dir, and --freesurfer_dir are all required. Each of these should be set to the overall directory (i.e. not an individual subject's directory). 

    sudo docker run -it --rm \
    -v /data/HCP_Raw:/data/HCP_Raw/ \
    -p 8888:8888 \
    winstonian3/connectome_harmonic:latest \
    /data/HCP_Raw/derivatives \
    participant \
    --mrtrix_dir /data/HCP_Raw/derivatives/MRtrix3_connectome-preproc \
    --ciftify_dir /data/HCP_Raw/derivatives/ciftify \
    --freesurfer_dir /data/HCP_Raw/derivatives/freesurfer \
    --participant_label 105923

Below are the full list of options and their descriptions. Please note that CHAP is a rough draft and none of these options are guaranteed to work at this time. Please report any errors or issues on Github, and thank you for your patience.

        usage: entrypoint_script.py [-h]
                                    [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                                    [--mrtrix_dir MRTRIX_DIR] [--hcp_dir HCP_DIR]
                                    [--ciftify_dir CIFTIFY_DIR]
                                    [--freesurfer_dir FREESURFER_DIR] [--evecs EVECS]
                                    [--nnum NNUM] [--tol TOL] [--skip_func SKIP_FUNC]
                                    [--diff_pipeline DIFF_PIPELINE]
                                    [--streamlines STREAMLINES]
                                    [--mask_med_wall MASK_MED_WALL]
                                    [--binarize BINARIZE]
                                    [--calculate_criticality CALCULATE_CRITICALITY]
                                    output_dir analysis_level

        Connectome Harmonic Analysis Pipeline (CHAP)

        positional arguments:
          output_dir            CHAP output directory (path)
          analysis_level        Participant or group mode. Only participant mode
                                supported for now.

        optional arguments:
          -h, --help            show this help message and exit
          --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                                Participant label(s) (not including sub-). If this
                                parameter is not provided all subjects will be
                                analyzed. Multiple participants can be specified with
                                a space separated list
          --mrtrix_dir MRTRIX_DIR
                                bids/mrtrix_connectome preproc output directory.
                                Required for CHAP-BIDS pipeline
          --hcp_dir HCP_DIR     HCP (min) preprocessed data directory. First level
                                should be test and retest folders OR if one session
                                just downloads. If test-retest, downloads go in
                                respective session folders. Required for CHAP-HCP
                                pipeline.
          --ciftify_dir CIFTIFY_DIR
                                Ciftify dir (required for CHAP-BIDS). Specify the
                                whole directory (i.e. not individual subject's
          --freesurfer_dir FREESURFER_DIR
                                Freesurfer dir (required for CHAP-BIDS). Specify the
                                whole directory (i.e. not individual subject's
          --evecs EVECS         Number of eigenvectors (harmonics) to compute. Default
                                is 100 (minus first trivial harmonic)
          --nnum NNUM           Number of nearest neighboring surface vertices to
                                assign to each streamline endpoint. Default = 60
          --tol TOL             (Tolerance) search radius of nearest neighbor search
                                for matching endpoints to surface vertices in mm.
                                Default = 1
          --skip_func SKIP_FUNC
                                Just find structural harmonics, no spectra.
          --diff_pipeline DIFF_PIPELINE
                                Choices: msmt_5tt pipeline or dhollander pipeline
                                based on bids/mrtrix3_connectome. Choose msmt or
                                dholl.
          --streamlines STREAMLINES
                                Number of streamlines in MRtrix tckgen
          --mask_med_wall MASK_MED_WALL
                                Mask out medial wall vertices. Default is True.
          --binarize BINARIZE   Binarize structural connectivity matrix. Default is
                                True
          --calculate_criticality CALCULATE_CRITICALITY
                                compute the criticality of the spectra across subjects

.. include:: links.rst

.. _Usage :

Usage Notes
===========

**Notes on fMRIprep**:

`fMRIprep`_ should be run first. Make
sure the "anat" output space is selected, and make sure FreeSurfer is
run.

**Notes on ciftify**:

`Ciftify`_ should be run with
the following two options: 1) --read-from-derivatives should point to
the derivatives directory that contains your fMRIprep and FreeSurfer
output folders. 2)--resample-to-T1w32k

**Notes on MRtrix3\_connectome**:

`MRtrix3_connectome`_
should be run with the following options: 1) the "preproc" level should
be run as opposed to the "participant" level. 2) --t1w\_preproc should
point to the location of the desc-preproc T1w in the fMRIprep output
directory.

Command-Line Arguments
----------------------
::

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
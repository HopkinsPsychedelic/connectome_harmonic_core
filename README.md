# CHAP (Connectome Harmonic Analysis Pipeline)
<img width="174" alt="Screen Shot 2021-08-24 at 9 29 19 PM" src="https://user-images.githubusercontent.com/61159065/130711344-9a354697-8525-4fea-a05a-4659f759e962.png"> <img width="174" alt="Screen Shot 2021-08-24 at 9 31 09 PM" src="https://user-images.githubusercontent.com/61159065/130711463-db54cb33-889b-4ecb-af33-b75cfb9d2930.png"> <img width="174" alt="Screen Shot 2021-08-24 at 9 30 38 PM" src="https://user-images.githubusercontent.com/61159065/130711430-df5b0c1f-38a8-4d43-8676-be640892c898.png"> <img width="174" alt="Screen Shot 2021-08-24 at 9 31 28 PM" src="https://user-images.githubusercontent.com/61159065/130711497-00af5938-36f1-4543-9f8f-eea8c1d3475f.png">

# Intro 
An outstanding challenge in neuroimaging research is interpreting functional brain data (collected using functional magnetic resonance imaging, or fMRI) through the bounds of structural connectivity (often estimated using diffusion-weighted and/or T1-weighted MRI.) Attempts have typically included separate pipelines of analysis for functional and structural data, with a final step of qualitative comparison between the spatial distribution of functional and structural connectivity, or some post-hoc quantitative comparison of the two. Other attempts at understanding function through structure include the use of a structure to generate a graph through which function is then understood (Bullmore & Sporns 2009). 

Many recent attempts to understand and interpret function rely on a dimensionality reduction step, involving parcellation of brain space that reduces the unit of analysis from voxel to region or network. Since certain sets of voxels have shown a high degree of temporal correlation under certain task and task-free conditions, grouping these voxels together into parcels potentially yields more functionally meaningful units of analysis than individual voxels. Voxel-wise data is typically projected through the set of parcels, at which point the temporal dynamics and/or connectivity profiles may be analyzed. 

Common dimensionality reduction techniques include atlas-based and independent component analysis (ICA) methods, which leverage either a-priori or data-driven functional relationships to interpret functional data. A limitation of these approaches, however, is that they do not explicitly account for structure-function relationships intrinsic to the structural connectivity of the brain. Similarly, correlative or graph-based methods of relating structural and functional data have historically come with a-priori assumptions that are built into the particular basis set that is used. 

Connectome harmonic decomposition, first described in Atasoy 2016, is an emerging technique for interpreting functional brain data through a basis set derived exclusively from brain structure (Atasoy, 2016). Connectome harmonics are defined as the eigenvectors of the graph laplacian computed on the structural connectome. Laplacian eigenfunctions describe a diversity of natural processes, including electron wave functions, standing waves formed by the vibrations of musical instruments, and even animal coat patterns. In the case of the connectome, the graph Fourier basis manifests as a set of patterns on the cortical surface—from now on referred to as connectome harmonics. Within this framework, functional data at each timepoint can be interpreted as a weighted sum of harmonics. Therefore, unlike atlas-based or ICA methods, this method uses structurally-derived parcels to describe functional data. 
 
In a study of 10 individuals, Atasoy and colleagues demonstrated that relatively few harmonics can reconstruct most canonical functional resting state networks as defined by the Yeo atlas. In particular, the Default Mode Network reliably showed a large degree of mutual information with the 9th harmonic (+/-2 due to individual differences). Generally, lower-order sensory and emotional processing networks—such as visual, somato-motor, and limbic—correlated significantly with low-frequency harmonics, whereas networks implicated in higher cognitive functions—such as control, dorsal attention, ventral attention, and default mode—recruited a wider spectrum of harmonics, including higher-frequency ones. 

One major line of research has focused on the spectral power—estimated via fMRI—of high versus low frequency harmonics under different levels of conscious arousal. In 2017, Atasoy et. al demonstrated that ingestion of the 5-HT2A agonist LSD increased the relative power of high frequency harmonics compared to low-frequency harmonics in a placebo-controlled study. Under anesthesia, Luppi and colleagues showed the opposite effect—a shift in power from high-frequency to low-frequency harmonics. This pattern was also observed in patients with disorders of consciousness, but only unconscious (fMRI-) patients rather than covertly conscious (fMRI+) patients, further suggesting that higher-frequency harmonics may be functionally relevant to high-level cognition. 

# CHAP Overview
Connectome harmonic decomposition is a relatively new framework; in addition to theoretical considerations of the functional significance of harmonic modes, there are myriad basic, methodological questions that require answering. Although a growing number of groups are employing the connectome harmonic paradigm to interpret fMRI, and more recently EEG, data, there currently exists no open-source software program to perform this analysis. 

To fill this need, we created a containerized and open-source software pipeline, CHAP (Connectome Harmonic Analysis Pipeline). Using functions from several top open-source neuroimaging libraries such as FreeSurfer, MRtrix3, FSL, and fMRIprep, CHAP can compute harmonics on any dataset containing adequate structural and diffusion data. CHAP also projects BOLD fMRI data through harmonics to output a variety of functional spectra that are useful in analyzing resting state and event-related timeseries data.

CHAP contains two distinct arms: CHAP-HCP and CHAP-BIDS, which refer to the type of data that each accepts. 

CHAP-HCP takes as input Human Connectome Project data preprocessed according to the HCP minimal preprocessing pipeline. For each set of data from each subject and session, data contained in the HCP Structural, Diffusion, and Structural Extended downloads are used to compute harmonics. HCP resting state and task-based fMRI data are optional; when present, CHAP-HCP will output functional spectra for these scans. 

CHAP-BIDS can compute connectome harmonics on any dataset containing T1w and diffusion data with reverse phase-encoding fieldmaps. At this time, this arm requires the outputs of three open-source BIDS (Brain Imaging Data Structure) apps: bids/mrtrix3_connectome, for preprocessing of diffusion data, nipreps/fmriprep, for preprocessing of fMRI data and running FreeSurfer, and tigrlab/fmriprep_ciftify, for resampling of surfaces and functional timeseries. 

More details on how CHAP works can be found at the bottom of this page. 

# Instructions on usage

**Notes on fMRIprep**:

[fMRIprep](https://fmriprep.org/en/stable/) should be run first. Make sure the "anat" output space is selected and make sure FreeSurfer is run.

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

In the below example, I run CHAP-BIDS. In addition to the two positional arguments discussed above, --mrtrix_dir, --ciftify_dir, and --freesurfer_dir are all required. Each of these should be set to the overall output directory of that software (i.e. not an individual subject's directory). 

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

Below is the full list of options and their descriptions. Please note that CHAP is a rough draft and none of these options are guaranteed to work at this time. Please report any errors or issues on Github, and thank you for your patience.

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

# More on CHAP

**Connectomes Overview**:

In CHAP, we represent the structural connectome with a 59,412x59,412 adjacency matrix, an unprecedented level of resolution for this paradigm. We derive local connections from physical adjacency on the cortical surface and estimate long range connections between vertices using diffusion tractography.

**Surface Matrices**:

CHAP uses white matter surfaces resampled to 32k native space. In CHAP-HCP, preprocessed surface GIfTI files are retrieved from the HCP Structural download. In CHAP-BIDS, fmriprep is used to run FreeSurfer. Ciftify then resamples native FreeSurfer surfaces to 32k space. In both pipelines, surface files are located at:
{subject}/T1w/fsaverage_LR32k/{subject}.{hem}.white.32k_fs_LR.surf.gii Vertex coordinate and connectivity information from surface GIfTIs is loaded using functions from the Nibabel Python library. The left and right hemispheres are concatenated to generate an unmasked sparse adjacency matrix of size 64,948x64,948 in which connections between physically contiguous vertices are represented by a 1 and non-connections by a 0. 

**Diffusion preprocessing**:

CHAP HCP:

In CHAP-HCP, all processing of diffusion data subsequent to the HCP minimal preprocessing pipeline is performed using tools contained in the open-source software package MRtrix3. First, a tissue segmentation for use in anatomically constrained response estimation, fiber orientation distribution (FOD) computation, and tractography is derived using MRtrix’s hybrid surface-volume segmentation (hsvs) algorithm. The hsvs algorithm relies on FreeSurfer output, which is retrieved from the HCP Structural Extended download. 
Next, a FOD is calculated from the preprocessed diffusion image, and an estimated response function is derived using the MRtrix3 dhollander algorithm. As a default, 10 million streamlines are generated from the FOD using the Second-order Integration over Fiber Orientation Distributions (iFOD2) probabilistic tractography algorithm with 250 mm as the maximum streamline length and the “power” parameter set to 0.33. The seeding and termination of tracks is constrained to the white matter/gray matter interface as defined by FreeSurfer segmentation using the “anatomically constrained tractography” option. The spatial coordinates of the endpoints of each streamline are then extracted.

CHAP-BIDS:

In CHAP-BIDS, we recommend users run the open-source BIDS app bids/MRtrix3_connectome for preprocessing of diffusion data. MRtrix3_connectome first performs bias field correction and brain masking of the T1w image. To preprocess the diffusion images, it does Gibbs ringing removal; motion, eddy current, and EPI distortion correction; outlier detection and replacement; brain masking, bias field correction, and intensity normalization, and finally, rigid-body registration and transformation to the T1w image. CHAP-BIDS accepts the outputs of MRtrix3_connectome and runs an identical tractography reconstruction pipeline to CHAP-HCP.

**Long-Range Connectivity Matrices (same between CHAP-BIDS and CHAP-HCP)**:

In order to encode the information contained in each subject’s tractogram into the surface matrix, CHAP generates a long-range connectivity matrix by associating each surface vertex with its 45 nearest streamline endpoints within 3 mm (configurable parameters). The other endpoint of each of these streamlines is associated with its nearest neighboring surface vertex. Nearest-neighbor computation is conducted using SciKit-Learn’s kd-tree nearest neighbors algorithm. This connectivity information is stored in an adjacency matrix of size 64,948x64,948 in which long-range connections between surface vertices are represented by a 1 and non-connections by a 0.

**Generating Structural Harmonics**:

We concatenate the cortical surface and long-range connectivity matrices to estimate a structural connectome graph with dimensions 64,948x64,948. We then use the HCP mask to mask out vertices located on the medial wall where streamline termination is anatomically implausible. Masking reduces the matrix to size 59,412x59,412. 

Using functions from the Scipy.sparse library, we then compute the eigenvectors and eigenvalues of the graph Laplacian of the connectome. By default, CHAP saves 100 eigenvectors (harmonics) in a Numpy array of size 59,412x100, where the first (zero-th) harmonic is the trivial solution. Visualization toolkit (.vtk) files containing the set of harmonics projected on the cortical surface are also saved. 

**More incoming...**

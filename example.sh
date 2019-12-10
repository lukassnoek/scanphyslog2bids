#!/usr/bin/env bash
scanphyslog2bids \
    --file ../sub-02/func/sub-02_task-face_run-1_physio.log \
    --fmri ../sub-02/func/sub-02_task-face_run-1_bold.nii.gz \
    --sf 496  # use 500 when using a wireless physio recorder
    --triggermethod gradient_log \
    --outdir ../sub-02/func \  # default
    --plottraces \
    --plotalignment \
    --derivdir ../derivatives/physiology  # for plots

### Note 
# You don't have to supply the associated fMRI file, but in that case
# you need to specify the `tr`, and `ndyns` parameters (which can be
# ignored if you supply the `fmri` parameter).


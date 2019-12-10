""" Example usage of Python interface of scanphyslog2bids.
This shows how to set up a parallelized conversion workflow. """

import numpy as np
import os.path as op
import nibabel as nib
from glob import glob
from scanphyslog2bids.core import CouldNotFindThresholdError, PhilipsPhysioLog
from joblib import Parallel, delayed


def _run_parallel(log):
    """ Function for Parallel call """
    TRIGGER_METHOD = 'gradient_log'

    nii = log.replace('_physio.log', '_bold.nii.gz')
    vols = nib.load(nii).shape[-1]
    tr = np.round(nib.load(nii).header['pixdim'][4], 3)
    print(f'\nProcessing {log}: dyns={vols}, TR={tr:.3f}, method={TRIGGER_METHOD}')
    
    try:
        phlog = PhilipsPhysioLog(f=log, tr=tr, n_dyns=vols, sf=496, manually_stopped=False)  # init
        phlog.load()
        phlog.align(trigger_method=TRIGGER_METHOD)  # load and find vol triggers
        out_dir = op.join(f'../derivatives/physiology/figures')
        phlog.plot_alignment(out_dir=out_dir)  # plots alignment with gradient
        phlog.to_bids()  # writes out .tsv.gz and .json files
        phlog.plot_traces(out_dir=out_dir)
        phlog.plot_alignment(out_dir=out_dir)
    except CouldNotFindThresholdError:  # something went wrong with gradient-based alignment
        print(f"Could not find threshold for {log}")


logs = sorted(glob('../sub-*/func/*_physio.log'))
Parallel(n_jobs=10)(delayed(_run_parallel)(log) for log in logs)

import joblib as jl
import nibabel as nib
from glob import glob
from main import CouldNotFindThreshold, PhilipsPhysioLog

logs = sorted(glob('sub-*/func/*SCANPHYSLOG*.log'))

def _run_parallel(log):
    sub_name = op.basename(log).split("_")[0]
    trigger_method = 'gradient_log'  # or: 'interpolation', 'vol_triggers'
    nii = log.replace('_recording-respcardiac_physio.txt', '_bold.nii.gz')
    vols = nib.load(nii).shape[-1]
    tr = np.round(nib.load(nii).header['pixdim'][4], 3)
    print(f'\nProcessing {log}: dyns={vols}, TR={tr:.3f}, method={trigger_method}')
    
    try:
        phlog = PhilipsPhysioLog(f=log, tr=tr, n_dyns=vols, sf=496, manually_stopped=False)  # init
        phlog.load().align(trigger_method=trigger_method)  # load and find vol triggers
        out_dir = op.join(f'derivatives/physiology/{sub_name}')
        phlog.plot_alignment(out_dir=out_dir)  # plots alignment with gradient
        phlog.to_bids()  # writes out .tsv.gz and .json files
        to_return = None
    except CouldNotFindThresholdError:  # something went wrong with gradient-based alignment
        print(f"Could not find threshold for {log}")
        to_return = log

    return to_return

    wrong = jl.Parallel(n_jobs=10)(jl.delayed(_run_parallel)(log) for log in logs)
    wrong = [w for w in wrong if w is not None]  # filter wrong conversions
    np.savetxt('derivatives/physiology/wrong_conversions.txt', wrong, fmt='%s')
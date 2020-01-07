import os
import shutil
import numpy as np
import os.path as op
import nibabel as nib
from glob import glob
from scanphyslog2bids.core import PhilipsPhysioLog


def test_python_interface():
    here = op.dirname(__file__)
    data_dir = op.join(here, 'data')
    logs = [
        (op.join(data_dir, 'example_for_gradient_log.log'), 'gradient_log', 496),
        (op.join(data_dir, 'example_for_interpolation.log'), 'interpolate', 496),
        (op.join(data_dir, 'example_for_vol_markers.log'), 'vol_markers', 500)
    ]

    output_dir = op.join(data_dir, 'derivatives', 'physiology')
    os.makedirs(output_dir, exist_ok=True)


    for (f, method, sf) in logs:

        fmri_img = nib.load(f.replace('.log', '.nii.gz'))
        ndyns = fmri_img.shape[-1]
        tr = np.round(fmri_img.header['pixdim'][4], 3)

        print(f'\nProcessing {f}: sf={sf}, dyns={ndyns}, TR={tr:.3f}, method={method}')
        phlog = PhilipsPhysioLog(f=f, tr=tr, n_dyns=ndyns, sf=sf, manually_stopped=False)
        phlog.load()
        phlog.align(trigger_method=method)  # load and find vol triggers
        
        phlog.to_bids()  # writes out .tsv.gz and .json files
        phlog.plot_alignment(out_dir=output_dir)  # plots alignment with gradient
        phlog.plot_traces(out_dir=output_dir)

    # Teardown test
    if op.isdir(op.dirname(output_dir)):
        shutil.rmtree(op.dirname(output_dir))

    for ext in ['tsv.gz', '.json']:
        files = glob(op.join(data_dir, f'*{ext}'))
        _ = [os.remove(f) for f in files]

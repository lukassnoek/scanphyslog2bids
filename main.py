import os
import os.path as op
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from glob import glob
from copy import copy
import nibabel as nib
import gzip
import json


class CouldNotFindThresholdError(Exception):
    pass


class PhilipsPhysioLog:
    """ Reads, converts, and aligns Philips physiology files (SCANPHYSLOG).

    Work in progress!
    """
    def __init__(self, f, fmri_file=None, n_dyns=100, sf=496, tr=None, manually_stopped=False):
        """ Initializes PhilipsPhysioLog object. """
        self.f = f
        self.n_dyns = n_dyns
        self.sf = sf  # sampling freq
        self.tr = tr  # TR in secs
        self.manually_stopped = manually_stopped
        self.fmri_file = fmri_file

        if fmri_file is not None:
            img = nib.load(fmri_file)
            self.tr = img.header['pixdim'][4]
            self.n_dyns = img.header['dim'][4]

        if self.tr is None:
            raise ValueError("Please provide a TR")

        self.n_trig = n_dyns + 1 if manually_stopped else n_dyns
        self.trs = self.tr * self.sf

    def load(self):

        with open(self.f, 'r') as f_in:
            for i, line in enumerate(f_in):
                if line[:2] != '##':
                    break

            txt = f_in.readlines()
            txt = [line.replace('  ', ' ').replace('\n', '') for line in txt if line != '#\n']
            self.markers = np.array([s.split(' ')[9] for s in txt])

        m_start_idx = np.where(self.markers == '0100')[0]
        if len(m_start_idx) == 0:
            m_start_idx = 0
        else:
            m_start_idx = m_start_idx[-1]

        m_end_idx = np.where(self.markers == '0020')[0]
        if len(m_end_idx) == 0:
            print("WARNING: No end marker in phys-file!")
            m_end_idx = len(txt) - 1
        else:
            m_end_idx = m_end_idx[-1]

        self.m_start_idx = m_start_idx
        self.m_end_idx = m_end_idx
        self.dat = np.loadtxt(self.f, dtype=int, usecols=np.arange(9))
        self.n = self.dat.shape[0]
        self.grad = self.dat[:, (6, 7, 8)]

        return self

    def align(self, trigger_method='gradient_log', which_grad='y'):

        found_start_of_grad = False
        custom_end_idx = copy(self.m_end_idx)
        
        while not found_start_of_grad:

            if self.grad[custom_end_idx, :].any():
                found_start_of_grad = True
            else:
                custom_end_idx -= 1

            if custom_end_idx < 0:

                if trigger_method == 'gradient_log':
                    print("WARNING: gradients were not logged, using interpolation")
                    trigger_method = 'interpolate'

                custom_end_idx = self.m_end_idx
                break

        self.c_end_idx = custom_end_idx

        if trigger_method == 'gradient_log':
            self._determine_triggers_by_gradient(which_grad)
        elif trigger_method == 'interpolate':
            self._determine_triggers_by_interpolation()
        elif trigger_method == 'vol_triggers':
            self._determine_triggers_by_volume_markers()

        self.trigger_diffs = np.diff(np.r_[self.real_triggers, self.c_end_idx])
        diff_last_vol = (self.c_end_idx - self.real_triggers[-1])
        if np.abs(diff_last_vol - self.trs) > 10:
            diff_in_sec = diff_last_vol / self.sf
            print(f"WARNING: last trigger has a duration of {diff_in_sec:.2f} seconds!")
        
        # Weird triggers = those with a diff much larger/smaller than TR
        # (except the last one)
        weird_triggers_idx = np.abs(self.trigger_diffs - self.trs) > 3
        if np.abs(diff_last_vol - self.trs) > 3:
            # The last vol should not counts as a weird trigger
            weird_triggers_idx[-1] = False

        self.weird_triggers = self.real_triggers[weird_triggers_idx]

        if self.weird_triggers.size > 0:
            weird_triggers_diff = self.trigger_diffs[weird_triggers_idx]        
            print(f"WARNING: found {self.weird_triggers.size} weird triggers with the following durations:")
            print(weird_triggers_diff)

        if self.manually_stopped:
            self.real_triggers = self.real_triggers[:-1]
        
        m_diff, std_diff = self.trigger_diffs[:-1].mean() / self.sf, self.trigger_diffs[:-1].std() / self.sf 
        print(f"INFO: Found {self.real_triggers.size} triggers with a mean duration of {m_diff:.5f} ({std_diff:.5f})!")

        return self

    def _check_for_extra_triggers(self, real_triggers):
        """ Removes 'extra' (erroneous) triggers """
        trigger_diffs = np.r_[np.diff(real_triggers), self.trs]
        extra_idx = np.abs(self.trs - trigger_diffs) > 5
        prob_extra = np.diff(np.r_[extra_idx.astype(int), 0]) == -1
        if prob_extra.sum() > 0:
            for idx in np.where(prob_extra)[0]:
                pre, post = trigger_diffs[idx-1], trigger_diffs[idx]
                if np.abs((pre + post) - (self.trs * 2)) < 10:
                    print(f"WARNING: Interpolated extra trigger: {(pre / self.trs):.2f} (pre) and {(post / self.trs):.2f} (post)")
                    real_triggers = real_triggers[real_triggers != real_triggers[idx]]
                    real_triggers = np.r_[real_triggers, real_triggers[idx-1] + self.trs]

            real_triggers.sort()
        return real_triggers

    def _check_for_missed_triggers(self, real_triggers):

        trigger_diffs = np.diff(np.r_[real_triggers, self.c_end_idx])
        prob_missed = real_triggers[np.abs(trigger_diffs - (self.trs * 2)) < 10]
        for i, trig in enumerate(prob_missed):
            real_triggers = np.r_[real_triggers, trig + self.trs]

        real_triggers.sort()
        return real_triggers

    def _determine_triggers_by_volume_markers(self):
        # remove markers before start and after end
        init_triggers = np.where(np.logical_or(self.markers == '0200', self.markers == '0202'))[0]
        init_triggers = self._check_for_missed_triggers(init_triggers)
        init_triggers = self._check_for_extra_triggers(init_triggers)

        real_triggers = init_triggers[-self.n_trig:]
        if len(real_triggers) != self.n_trig:
            raise ValueError(f"ERROR: expected to find {self.n_dyns} triggers, but found {len(real_triggers)}"
                             f"(and {len(init_triggers)} init triggers)")

        self.real_triggers = real_triggers

    def _determine_triggers_by_interpolation(self):
        
        if self.manually_stopped:
            print("WARNING: using the interpolation method with manually stopped scans is a very bad idea!")
        
        OFFSET = int((166./500) * self.sf)
        assumed_start = int(self.c_end_idx - OFFSET - (self.trs * self.n_trig))
        # determine samples per tr
        leftover = np.round(self.trs % 1, 1)
        p = int(1 / (1 - leftover))
        diffs = np.array([int(np.floor(self.trs)) if i % p == 0
                          else int(np.ceil(self.trs)) for i in range(self.n_dyns)])
        diffs[0] += assumed_start
        self.real_triggers = np.cumsum(diffs)

    def _determine_triggers_by_gradient(self, which_grad):

        self.align_grad = self.grad[:, {'x': 0, 'y': 1, 'z': 2}[which_grad]]
        # set prescan stuff to zero
        self.approx_start = self.c_end_idx - (self.sf * self.tr * self.n_trig) - self.trs * 0.05
        grad = self.align_grad.copy()
        grad[np.arange(self.n) < self.approx_start] = 0
        
        thr = self.align_grad.min()
        while True:
            # Find potential triggers
            real_triggers = np.where(grad < thr)[0].astype(int)
            trigger_diffs = np.diff(np.r_[real_triggers, self.c_end_idx])
            real_triggers = real_triggers[trigger_diffs > 2]  # remove doubles
            
            # Check for missed triggers
            if (self.n_dyns - real_triggers.size) < 10:
                real_triggers = self._check_for_missed_triggers(real_triggers)
                real_triggers = self._check_for_extra_triggers(real_triggers)

            if real_triggers.size == self.n_trig:
                break # Found it!

            thr += 1
            if thr > 0:
                raise CouldNotFindThresholdError("Could not find threshold!")

        self.real_triggers = real_triggers

    def to_bids(self):
        print((self.m_end_idx - self.real_triggers[-1] - self.trs))
        print((self.c_end_idx - self.real_triggers[-1] - self.trs))
        base_name, _ = op.splitext(self.f)

        time = np.arange(self.n) / self.sf
        start = time[self.real_triggers[0]]
        time = time - start

        info = {
           "SamplingFrequency": self.sf,
           "StartTime": time[0],
           "Columns": ["cardiac", "respiratory", "trigger"]
        }

        with open(f'{base_name}.json', "w") as write_file:
            json.dump(info, write_file, indent=4)

        data = self.dat[:, 4:6]
        pulses = np.zeros(self.n)
        pulses[self.real_triggers] = 1
        data = np.c_[data, pulses]
        tsv_out = f'{base_name}.tsv'
        np.savetxt(tsv_out, data, delimiter='\t')
        with open(tsv_out, 'rb') as f_in, gzip.open(tsv_out + '.gz', 'wb') as f_out:
            print(f"INFO: Saving to {tsv_out} ...")
            f_out.writelines(f_in)
        os.remove(tsv_out)

    def plot_alignment(self, win=4000, out_dir=None):

        n_weird = len(self.weird_triggers)
        if n_weird > 5:
            print(f"WARNING: found {n_weird} weird triggers! Only going to plot the first 5")
            n_weird = 5

        fig, ax = plt.subplots(nrows=4 + n_weird, figsize=(30, 9 + 3 * n_weird))

        if hasattr(self, 'align_grad'):
            amp = self.align_grad.min() * .25
        else:
            amp = 1

        real_triggers = np.zeros(self.n)
        self.real_triggers = self.real_triggers.astype(int)
        real_triggers[self.real_triggers] = amp
        start_end = np.zeros(self.n)
        start_end[self.m_start_idx] = amp
        start_end[self.m_end_idx] = amp
        start_end[self.c_end_idx] = amp * 2

        windows = [
            ('Full', np.arange(self.n)[self.m_start_idx:]),
            ('Start', np.arange(self.m_start_idx, self.m_start_idx + win)),
            ('End', np.arange(self.m_end_idx - win, self.m_end_idx)),
        ]

        for weird_trig in self.weird_triggers[:n_weird]:
            trig_nr = np.where(weird_trig == self.real_triggers)[0][0]
            windows.append(
                (f'Weird trigger (#{trig_nr+1}) at {weird_trig}',
                 np.arange(weird_trig - win, np.min([weird_trig + win, self.c_end_idx])).astype(int)
            ))

        for i, (title, period) in enumerate(windows):

            ext_space = 300 if i == 0 else 50
            lw_trigs = 1 if i == 0 else 2
            legend = []
            if hasattr(self, 'align_grad'):
                ax[i].plot(period, self.align_grad[period], lw=0.5)
                legend.append('grad')

            ax[i].plot(period, start_end[period], lw=2)
            legend.append('start/end')

            ax[i].plot(period, real_triggers[period], lw=lw_trigs)
            ax[i].set_title(title, fontsize=15)
            ax[i].set_xlim(period[0] - ext_space, period[-1] + ext_space)
            ax[i].legend(legend + ['triggers'])
            if i > 2:
                ax[i].axvline(self.weird_triggers[i-3], ls='--', c='k', lw=0.8)
                
            if hasattr(self, 'approx_start') and i == 0:
                ax[0].axvline(self.approx_start, ls='--', c='k', lw=0.5)

        ax[-1].plot(self.trigger_diffs[:-1])
        ax[-1].set_title('Number of samples between triggers', fontsize=15)
        m_diff = np.mean(self.trigger_diffs[:-1])
        std_diff = np.std(self.trigger_diffs[:-1])
        ax[-1].text(0, self.trigger_diffs[:-1].max(), f'm = {m_diff:.0f} ({(m_diff / self.sf):.3f}), std = {std_diff:.0f} ({(std_diff / self.sf):.3f})', fontsize=20)
        fig.tight_layout()

        if out_dir is not None:
            out_dir = op.join(out_dir, 'figures')
            if not op.isdir(out_dir):
                os.makedirs(out_dir, exist_ok=True)

            f_out = op.join(out_dir, op.splitext(op.basename(self.f))[0] + '_alignment.png')
            if n_weird > 0:
                f_out = f_out.replace('.png', '_CHECKOUT.png')  # for debugging purposes
            fig.savefig(f_out, dpi=100)
            plt.close()

if __name__ == '__main__':
    import nibabel as nib
    import joblib as jl
    from glob import glob
    logs = sorted(glob('../sub-*/func/*work*.txt'))

    def _run_parallel(log):
        sub_name = op.basename(log).split("_")[0]
        trigger_method = 'gradient_log'
        nii = log.replace('_recording-respcardiac_physio.txt', '_bold.nii.gz')
        vols = nib.load(nii).shape[-1]
        tr = np.round(nib.load(nii).header['pixdim'][4], 3)
        print(f'\nProcessing {log}: dyns={vols}, TR={tr:.3f}, method={trigger_method}')
        try:
            ms = True if 'stopsignal' in log else False
            l = PhilipsPhysioLog(f=log, tr=tr, n_dyns=vols, sf=496, manually_stopped=ms)
            l.load().align(trigger_method=trigger_method)
            out_dir = op.join(f'../derivatives/physiology/{sub_name}')
            l.plot_alignment(out_dir=out_dir)
            l.to_bids()
            to_return = None
        except CouldNotFindThresholdError:
            print(f"Could not find threshold for {log}")
            to_return = log

        return to_return

    wrong = jl.Parallel(n_jobs=10)(jl.delayed(_run_parallel)(log) for log in logs)
    wrong = [w for w in wrong if w is not None]
    #for log in logs:
    #    _run_parallel(log)
    np.savetxt('../derivatives/physiology/wrong_conversions.txt', wrong, fmt='%s')

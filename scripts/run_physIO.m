function [  ] = run_retroicor(phys_file, fmri_file, sf, save_dir, per_slice)
% Runs TAPAS PhysIO RETROICOR/HRV/RVT estimation

% load sf from json and tr/n_slices/n_vols from nii
% assumes that the following "package" is installed:
% https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
nii = load_untouch_header_only(fmri_file);
n_slices = nii.dime.dim(2);
n_vols = nii.dime.dim(5);
tr = nii.dime.pixdim(5);

physio = tapas_physio_new();

[~,basename,~] = fileparts(phys_file);
ricor_out = strrep(basename,'physio.tsv','desc-retroicor_regressors.tsv');
mat_out = strrep(basename,'physio.tsv','desc-PhysIO_info.mat');

fprintf('Trying to create %s with %i vols, %i slices, and TR=%.3f...', ricor_out, n_vols, n_slices, tr);

% Turn off plots
set(0,'DefaultFigureVisible','off');

%% Individual Parameter settings. Modify to your need and remove default settings
physio.save_dir = {save_dir};
physio.log_files.vendor = 'BIDS';
physio.log_files.cardiac = {phys_file};
physio.log_files.align_scan = 'first';
physio.log_files.sampling_interval = 1/sf;
physio.scan_timing.sqpar.Nslices = n_slices;
physio.scan_timing.sqpar.TR = tr;
physio.scan_timing.sqpar.Ndummies = 0;  % already removed
physio.scan_timing.sqpar.Nscans = n_vols;

if per_slice == 1
    physio.scan_timing.sqpar.onset_slice = 1:n_slices;
else
    physio.scan_timing.sqpar.onset_slice = floor(n_slices / 2);
end    

physio.scan_timing.sync.method = 'scan_timing_log';
physio.preproc.cardiac.modality = 'PPU';
physio.preproc.cardiac.initial_cpulse_select.method = 'auto_matched';
physio.preproc.cardiac.initial_cpulse_select.file = 'initial_cpulse_kRpeakfile.mat';
physio.preproc.cardiac.initial_cpulse_select.min = 0.4;
physio.preproc.cardiac.posthoc_cpulse_select.method = 'off';
physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80;
physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;
physio.model.orthogonalise = 'none';
physio.model.censor_unreliable_recording_intervals = false;
physio.model.output_multiple_regressors = ricor_out;
physio.model.output_physio = mat_out;
physio.model.retroicor.include = true;

physio.model.retroicor.order.c = 3;
physio.model.retroicor.order.r = 4;
physio.model.retroicor.order.cr = 1;
physio.model.rvt.include = true;
physio.model.hrv.include = true;    

physio.model.rvt.delays = 1;
physio.model.hrv.delays = 1;

physio.model.noise_rois.include = false;
physio.model.movement.include = false;
physio.model.other.include = false;
physio.verbose.level = 0;
physio.verbose.use_tabs = true;
physio.ons_secs.c_scaling = 1;
physio.ons_secs.r_scaling = 1;

% RETROICOR cardiac regressors [2 x 3 nOrderCardiac] = 6
% RETROICOR respiratory regressors [2 x 4 nOrderRespiratory] = 8
% RETROICOR cardXResp interaction regressors [4 x 1 nOrderCardiacXRespiratory] = 4
% HRV [nDelaysHRV] = 1
% RVT [nDelaysRVT] = 1 

%% Run physiological recording preprocessing and noise modeling
physio = tapas_physio_main_create_regressors(physio);

end


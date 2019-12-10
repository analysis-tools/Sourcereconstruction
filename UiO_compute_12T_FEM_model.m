function [headmodel_vol, headmodel_elec, lead_field_gray, lead_field] = UiO_compute_12T_FEM_model(vol, elec_aligned, save_folder, SegGray)

%% prepare sens (seems to be done in ft_prepare_leadfield anyway...)

% channels = ft_channelselection({'EEG'}, elec_aligned);
disp('Now doing bookkeeping to make sure volume and sensors are compatible (maybe?).')
disp('this will take many minutes. Probably hours. maybe around three of them? or less?')
disp('we think 3 hours is a good estimate')
[vol, elec_aligned] = ft_prepare_vol_sens(vol, elec_aligned, 'channel', elec_aligned.label);
%[volGray, elec_aligned] = ft_prepare_vol_sens(volGray, elec_aligned, 'channel', channels);

disp('saving models')
headmodel_vol = [save_folder '\headmodel_12T_FEM_prepared_sens_vol.mat'];
headmodel_elec = [save_folder '\headmodel_12T_FEM_prepared_sens_elec.mat'];

save(headmodel_vol,'vol','-v7.3');
save(headmodel_elec,'elec_aligned','-v7.3');
disp('headmodel and electrodes saved');


%% sourcemodel
disp('creating source model');
cfg                = [];
cfg.mri            = SegGray;
cfg.elec = elec_aligned;
cfg.vol = vol;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 6; %change resolution???
cfg.spmversion = 'spm12';
gridGray           = ft_prepare_sourcemodel(cfg);%, vol, elec_aligned);


%% leadfield

disp('calculating leadfield matrix');
disp('Will take at least 1 hour');
cfg = [];
cfg.spmversion = 'spm12';
cfg.elec = elec_aligned;
cfg.headmodel = vol;
cfg.grid = gridGray;
% cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
cfg.reducerank = 3; % or 3?
[grid,cfg] = ft_prepare_leadfield(cfg);

lead_field_gray = [save_folder '\leadfield_12T_FEM_gray-only.mat'];
lead_field = [save_folder '\LFcfg_12T_FEM.mat'];

save(lead_field_gray,'grid','-v7.3');
save(lead_field,'cfg','-v7.3');
disp('saved leadfield and logFile');
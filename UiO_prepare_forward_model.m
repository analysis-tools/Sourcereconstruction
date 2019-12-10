function [mesh_name,elec_name,vol_name,SegGray_name] = UiO_prepare_forward_model(mri_data, save_folder, loc_file, loc_path)
% for debugging
% mri_data = [subjID 'MIDA_withoutair.nii'];
% mri_path = save_folder;
% mri_data = [mri_path,mri_data];

%% load MIDA and segment
MIDA_mri = ft_read_mri(mri_data);

% cfg = [];
% cfg.spmversion     = 'spm12';
% cfg.resolution = 1;
% cfg.dim = [256 256 256];
% MIDA = ft_volumereslice(cfg,MIDA);

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[MIDA_mri] = ft_volumerealign(cfg,MIDA_mri);

% Define tissue types
Tissues = {'skin','eyes','muscle','fat','spongybone','compactbone','gray','cerebellargray','white','cerebellarwhite','csf','brainstem'};

skin = [1,33:35,37,39,51,85,86];
eyes = [55:59];
muscle = [38,42,60,61,63:84,88:96,98];
fat = [43,62];
spongybone = [52];
compactbone = [36,40,41,44:49,53,54,87];
gray = [3:5,7,8,10,16,17,19,20,21,99,116];
cerebellargray = [2];
white = [12,18,22,23,100:115];
cerebellarwhite = [9];
csf = [6,24,25,32];
brainstem = [11,13:15];

% special tissue for GrayMatter
SegGrayMatter = sort([gray cerebellargray]);

bg = 50;
bg = ismember(MIDA_mri.anatomy,bg);
MIDA_mri.anatomy(bg) = 0;


%% use warped MIDA to define segmented tissues

for i = 1:length(Tissues)
    Segment = ismember(MIDA_mri.anatomy,eval(Tissues{i}));
    eval([Tissues{i} ' = Segment;']);
end

Seg = struct;
Seg.skin = skin;
Seg.eyes = eyes;
Seg.muscle = muscle;
Seg.fat = fat;
Seg.spongybone = spongybone;
Seg.compactbone = compactbone;
Seg.gray = gray;
Seg.cerebellargray = cerebellargray;
Seg.white = white;
Seg.cerebellarwhite = cerebellarwhite;
Seg.csf = csf;
Seg.brainstem = brainstem;

Seg.dim = MIDA_mri.dim;
Seg.coordsys = MIDA_mri.coordsys;
Seg.transform = MIDA_mri.transform;
Seg.unit = MIDA_mri.unit;
Seg.transformorig = MIDA_mri.transformorig;

%% Seg Gray matter

grayMatter = ismember(MIDA_mri.anatomy,SegGrayMatter);

SegGray = struct;
SegGray.gray = grayMatter;
SegGray.dim = MIDA_mri.dim;
SegGray.coordsys = MIDA_mri.coordsys;
SegGray.transform = MIDA_mri.transform;
SegGray.unit = MIDA_mri.unit;
SegGray.transformorig = MIDA_mri.transformorig;


%% compute hexahedral meshes
cfg = [];
cfg.spmversion = 'spm12';
cfg.tissues = Tissues;
cfg.shift  = 0;
cfg.method = 'hexahedral';
cfg.downsample = 1;
%cfg.smooth = 'no';
mesh = ft_prepare_mesh(cfg,Seg);


%% headmodel
cfg        = [];
cfg.method ='simbio';
cfg.conductivity = [0.4348 0.5 0.1 0.04 0.04 0.0063 0.3333 0.2564 0.1429 0.1099 1.5385 0.1538];   % order follows mesh.tissuelabel
vol = ft_prepare_headmodel(cfg, mesh);   

%% Read sensor locations (Currently somewhat of a hack)
path_loc_spec = [loc_path,loc_file];
%path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];
%path_loc_attload = path_loc_spec; %Path to file
elec_xyz = ft_read_sens(path_loc_spec);


%% Adjust final position of electrodes by eye
% this is not necessary when using Simbio since it takes the closest vertex
% of the outer skin as electrode positions...
% 
disp('preparing for final adjustments of the electrode positions')
disp('this takes a couple of minutes')
cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_xyz;
cfg.headshape = vol;
elec_aligned  = ft_electroderealign(cfg);
elec_aligned = ft_convert_units(elec_aligned,'mm');

mesh_name = [save_folder '\mesh.mat'];
elec_name = [save_folder '\elec_aligned.mat'];
vol_name = [save_folder '\vol.mat'];
SegGray_name = [save_folder '\SegGray.mat'];

save(mesh_name,'mesh','-v7.3');
save(elec_name,'elec_aligned','-v7.3');
save(vol_name,'vol','-v7.3');
save(SegGray_name,'SegGray','-v7.3');

end


% %% Plot electrode positions after realignment 
% figure;
% hold on;
% ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
% % camlight
% % electrodes
% ft_plot_sens(elec_aligned,'label','label');
% title('Electrode positions after realignment');


% %%%%%% hard coded passage (flip channels for UiO eeg-cap; look for EM and
% %%%%%% EO for Sebastian eeg files)
% 
% %% old for UiO
% % mark and remove EOG channels
% idx_eog = cellfun(@(x)contains(x,'EO','IgnoreCase',true),elec_xyz.label) + cellfun(@(x)contains(x,'EM','IgnoreCase',true),elec_xyz.label);
% channels = [elec_xyz.label(~idx_eog)];
% 
% elec = ft_read_sens(path_loc_spec); %Read layout file
% 
% idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
% %Updating fields
% 
% elec.chanpos = elec.chanpos(idx_keep,:);
% elec.chantype = elec.chantype(idx_keep);
% elec.chanunit = elec.chanunit(idx_keep);
% elec.elecpos = elec.elecpos(idx_keep,:);
% elec.label = elec.label(idx_keep);
% 
% % realign channel order to subject channel order
% realignIdx = [];
% for i = 1:numel(elec.label)
%     realignIdx = [realignIdx, find(strcmp(upper(elec.label(:)),upper(channels(i))))];
% end
% 
% elec.label = elec.label(realignIdx);
% elec.chanpos = elec.chanpos(realignIdx,:);
% elec.chantype = elec.chantype(realignIdx);
% elec.chanunit = elec.chanunit(realignIdx);
% elec.elecpos = elec.elecpos(realignIdx,:);
% 
% logFile.EEGchanlabels = elec.label;
% 
% elec_aligned = elec;
% 
% %% new for Sebastian data
% 
% % idx_eog = cellfun(@(x)contains(x,{'EO', 'EM'},'IgnoreCase',true),elec_attload.label);
% % 
% % channels = [elec_attload.label(~idx_eog)];
% % 
% % elec = elec_attload;
% % 
% % 
% % 
% % idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
% % %Updating fields
% % 
% % elec.chanpos = elec.chanpos(idx_keep,:);
% % elec.chantype = elec.chantype(idx_keep);
% % elec.chanunit = elec.chanunit(idx_keep);
% % elec.elecpos = elec.elecpos(idx_keep,:);
% % elec.label = elec.label(idx_keep);
% % 
% % elec_aligned = elec;
% % elec_aligned = ft_convert_units(elec_aligned,'mm');
% % 
% % elec_aligned.elecpos = [elec_aligned.elecpos(:,2)*-1 elec_aligned.elecpos(:,1) elec_aligned.elecpos(:,3)];
% % elec_aligned.chanpos = [elec_aligned.chanpos(:,2)*-1 elec_aligned.chanpos(:,1) elec_aligned.chanpos(:,3)];
% 
% % elec = ft_determine_coordsys(elec_aligned);
% %%%%%%%

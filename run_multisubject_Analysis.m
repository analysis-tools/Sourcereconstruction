
%% Defining variables and parameters for source recontruction
% Creating variables
N = 2; % number of subjects
S = 1; % number of sessions
project_data = struct(); % structure for project data
participant_data = struct(); % structure for storing all data

participant_IDs = {'101','111','112','113','121','122','131','132','133',...
    '141','42','143','151','152','153','161','162','163','171','172','173',...
    '181','182','183','191','192'};
 
path_project = uigetdir([],'Choose a folder to store all data from this project!'); % should be in lagringshotell

while 1
    prompt = 'Do you want to use the standard folder structure? [y/n]: ';
    folder_structure = input(prompt,'s');
    
    if strcmp(folder_structure,'n') || strcmp(folder_structure,'N')
        
        % prompt user for relevant folder directories
        initalpath= uigetdir([],'Specify which folder you want the folder searches to start from');
        path_ft = uigetdir(initalpath,'Give me the Field Trip folder!'); % should be in lagringshotell
        path_spm = uigetdir(initalpath,'Give me the SPM12 folder!'); % should be in lagringshotell
        path_simbio = uigetdir(initalpath,'Give me the simbio folder (from fieldtrip)!'); % should be in lagringshotell
        [MIDA, MIDApath] = uigetfile('*.nii', 'Pick unmodified MIDA template'); % should be in lagringshotell
        [loc_file,loc_path] = uigetfile('*.*','Feed me the EEG electrode locations'); % should be in lagringshotell
        break
        
    elseif strcmp(folder_structure,'y') || strcmp(folder_structure,'Y')
        
        path_ft = 'Z:\source_reconstruction\fieldtrip-20191127\fieldtrip-20191127';
        path_spm = 'Z:\source_reconstruction\spm12\spm12';
        path_simbio = 'Z:\source_reconstruction\fieldtrip-20191127\fieldtrip-20191127\external\simbio';
        MIDA = 'MIDA_v1.nii';
        MIDApath = 'Z:\source_reconstruction\MIDA_model\MIDAv1.0\MIDA_v1.0\MIDA_v1_voxels\';
        loc_file = 'Oslo62_s.sfp';
        loc_path = 'Z:\source_reconstruction\';
        break
        
    else
        disp('answer only Y or N!')
    end
end

% storing relevant fields in data structure

project_data.path_project = [path_project,filesep];
project_data.MIDA = MIDA;
project_data.MIDApath = MIDApath;
project_data.loc_file = loc_file;
project_data.loc_path = loc_path;
project_data.path_ft = [path_ft, filesep];
project_data.path_spm = [path_spm, filesep];
project_data.path_simbio = [path_simbio, filesep];

%% Process MIDA model with fieldtrip ( I think this just needs to be done once )
restoredefaultpath;
addpath(path_ft);
ft_defaults;

while 1
    prompt = 'Is the MIDA model tissues already in your project folder? [y/n]: ';
    folder_structure = input(prompt,'s');
    
    if strcmp(folder_structure,'n') || strcmp(folder_structure,'N')
        
        tmp = matlab.desktop.editor.getActive;
        cd(fileparts(tmp.Filename));
        
        % processing the mida model for this project
        MIDAtissues = UiO_process_MIDA(MIDA, MIDApath, path_project);
        
        % storing relevant fields in data structure
        project_data.MIDAtissues = MIDAtissues;
        break
        
    elseif strcmp(folder_structure,'y') || strcmp(folder_structure,'Y')
        
        if exist([project_data.path_project 'MIDA_withoutair.nii']) == 2
            project_data.MIDAtissues.withoutair = [project_data.path_project 'MIDA_withoutair.nii'];
            project_data.MIDAtissues.gray = [project_data.path_project 'TPM_gray.nii'];
            project_data.MIDAtissues.white = [project_data.path_project 'TPM_white.nii'];
            project_data.MIDAtissues.soft = [project_data.path_project 'TPM_soft.nii'];
            project_data.MIDAtissues.bone = [project_data.path_project 'TPM_bone.nii'];
            project_data.MIDAtissues.CSF = [project_data.path_project 'TPM_CSF.nii'];
            project_data.MIDAtissues.background = [project_data.path_project 'TPM_background.nii'];
            break
        else
            folder_structure = 'n';
        end
        
    else
        disp('answer only Y or N!')
    end
end

disp('section completed')

%% Looping through all participants to prepare for FEM calculation
for n = 1:N

    
    subjID = participant_IDs{n};
    save_folder = [path_project, filesep, subjID, filesep]; % folder for participant
    if ~exist(save_folder, 'dir')
       mkdir(save_folder)
    end
    % storing relevant fields in data structure
    participant_data(n).subjID = subjID;
    participant_data(n).save_folder = save_folder;
    
    %% Creating (or loading) Nifti MRI of subject
    while 1
        prompt = ['Does a Nifti (.nii) version of subject ' subjID ' MRI exist [y/n]: '];
        nifti_exists = input(prompt,'s');
         
        if strcmp(nifti_exists,'n') || strcmp(nifti_exists,'N')
            restoredefaultpath;
            addpath(path_ft);
            ft_defaults;
            
            [subjectMRI_DICOM, subjectMRIpath_DICOM] = uigetfile([project_data.path_project '*.*'], 'Pick subject MRI in DICOM');
            
            disp('Loading MRI')
            mri = ft_read_mri([subjectMRIpath_DICOM subjectMRI_DICOM]);
            
            disp('Writing sources')
            mri_name = [subjID '_raw_MRI'];
            cfg = [];
            cfg.filename  = [save_folder filesep mri_name];
            cfg.filetype  = 'nifti';
            cfg.parameter = 'anatomy';
            ft_sourcewrite(cfg, mri);
            
            % [subjectMRI, subjectMRIpath] = uigetfile('*.nii', 'Pick subject MRI');
            subjectMRI = [mri_name '.nii'];
            subjectMRIpath = [save_folder filesep];
            disp('section completed')
            break
            
        elseif strcmp(nifti_exists,'y') || strcmp(nifti_exists,'Y')
            % When subject MRI available as NIFTI, start with this section
            [subjectMRI, subjectMRIpath] = uigetfile([project_data.path_project '*.nii'], 'Pick subject nifti file');
            disp('section completed')
            break
            
        else
            disp('answer only Y or N!')
        end
    end
    
    % storing relevant fields in data structure
    participant_data(n).subjectMRI = subjectMRI;
    participant_data(n).subjectMRIpath = subjectMRIpath;
end
% saving data structures
save([path_project, '\participant_data.mat'],'participant_data','-v7.3');
save([path_project, '\project_data.mat'],'project_data','-v7.3');

%% compute 12 tissue FEM model for each participant
if ~exist('participant_data','var')
    load([path_project, '\participant_data.mat']);
end
if ~exist('project_data','var')
    load([path_project, '\project_data.mat']);
end

%% Warp to subject space with SPM
for n = 1:N
    restoredefaultpath;
    addpath(genpath(path_spm));
    
    disp('This takes about 20 minutes. How about a coffee and a discussion about the nature of consciousness?')
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));
    
    subject_MIDA_tissues = UiO_transform_MIDA_to_subject(project_data.path_project, project_data.MIDAtissues.withoutair, participant_data(n).subjectMRI, participant_data(n).subjectMRIpath, subjID);
    % storing relevant fields in data structure
    participant_data(n).MIDA_tissues = subject_MIDA_tissues;
end
% saving data structures
save([path_project, '\participant_data.mat'],'participant_data','-v7.3');
save([path_project, '\project_data.mat'],'project_data','-v7.3');

%% compute 12 tissue FEM model for each participant
if ~exist('participant_data','var')
    load([path_project, '\participant_data.mat']);
end
if ~exist('project_data','var')
    load([path_project, '\project_data.mat']);
end

%% Prepare 12 tissue FEM model
disp('This takes about 10 minutes per subject and requires user input?')
for n = 1:N
    restoredefaultpath;
    addpath(path_ft);
    addpath(path_simbio);
    ft_defaults;
    
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));
    
    [mesh,elec_aligned,vol,SegGray] = UiO_prepare_forward_model(participant_data(n).MIDA_tissues.withoutair, save_folder, loc_file, loc_path);
    
    % storing relevant fields in data structure
    participant_data(n).mesh = mesh;
    participant_data(n).elec_aligned = elec_aligned;
    participant_data(n).vol = vol;
    participant_data(n).SegGray = SegGray;
    
end
% saving data structures
save([path_project, '\participant_data.mat'],'participant_data','-v7.3');
save([path_project, '\project_data.mat'],'project_data','-v7.3');

%% compute 12 tissue FEM model for each participant
if ~exist('participant_data','var')
    load([path_project, '\participant_data.mat']);
end
if ~exist('project_data','var')
    load([path_project, '\project_data.mat']);
end

disp('This takes about 8 hours per subject. Are you sure this is the best idea?')
for n = 1:N
    % loading data from 
    load(participant_data(n).vol);
    load(participant_data(n).elec_aligned);
    save_folder = participant_data(n).save_folder;
    load(participant_data(n).SegGray);

    [headmodel_vol, headmodel_elec, lead_field_gray, lead_field] = UiO_compute_12T_FEM_model(vol, elec_aligned, save_folder, SegGray);

    % storing relevant fields in data structure
    participant_data(n).headmodel_vol = headmodel_vol;
    participant_data(n).headmodel_elec = headmodel_elec;
    participant_data(n).lead_field_gray = lead_field_gray;
    participant_data(n).lead_field = lead_field;
    
end

clearvars vol elec_aligned SegGray

% saving data structures
save([path_project, '\participant_data.mat'],'participant_data','-v7.3');


%% Calculate inverse solution

% checking if data needs loading
if ~exist('participant_data','var')
    load([path_project, '\participant_data.mat']);
end
if ~exist('project_data','var')
    load([path_project, '\project_data.mat']);
end

% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));

% getting filenames for all EEG files (EEGlab structs)

prompt = 'Are your EEG files saved as EEGlab structs? [y/n]: ';
eegfile_prompt = input(prompt,'s');

if strcmp(eegfile_prompt,'y') || strcmp(eegfile_prompt,'Y')
    for n = 1:N
        eeg_data = struct([]);
        for s = 1:S
            [EEG_file,EEG_path] = uigetfile([project_data.path_project '*.*'],['Feed me the EEG data for subject (.mat eeglab struct)' participant_data(n).subjID ', session ' s]);
            eeg_data(s).file = EEG_file;
            eeg_data(s).path = EEG_path;
        end
        participant_data(n).eeg_data = eeg_data;
        
    end
    
else
    disp('Convert your files...');
end
% saving data structures
save([path_project, '\participant_data.mat'],'participant_data','-v7.3');

% checking if data needs loading
if ~exist('participant_data','var')
    load([path_project, '\participant_data.mat']);
end
if ~exist('project_data','var')
    load([path_project, '\project_data.mat']);
end

%
restoredefaultpath;
addpath(path_ft);
ft_defaults;

% calculating inverse solution for all data
for n = 1:N
    LF_data = participant_data(n).lead_field_gray;
    Head_data = participant_data(n).lead_field;
    save_folder = participant_data(n).save_folder;
    for s = 1:S
        
        MRI = ['r' participant_data(n).subjectMRI];
        MRIpath = participant_data(n).subjectMRIpath;
        EEG_file = participant_data(n).eeg_data(s).file;
        EEG_path = participant_data(n).eeg_data(s).path;
        ft = project_data.path_ft;
        UiO_inverse_model_eLORETA(MRI, MRIpath, EEG_file, EEG_path, LF_data, Head_data, save_folder, ft, save_folder);

    end
end


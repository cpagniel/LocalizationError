%% An experimental approach for reducing error in passive acoustic localizations 
%% of individual marine animals using in situ sound playbacks in complex 
%% shallow-water environments
% Pagniello 2021

% This code combines .mat file outputed from
% Pagniello_2021_localization_TDOA_**_replicas_SCOT.m. Change file name and
% output manually as needed.

%% Combine Files within Folders

rootdir = 'F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq9';
cd(rootdir);

files = dir('log*');
[~, reindex] = sort(str2double(regexp({files.name},'\d+','match','once')));
files = files(reindex);

range = 1:length(files);
files = files(range);

RL_rms_dB = NaN(length(files),1);
RL_pp_dB = NaN(length(files),1);

d_model = NaN(length(files),1);
d_combo = NaN(length(files),1);

SL_model_rms_dB = NaN(length(files),1);
SL_model_pp_dB = NaN(length(files),1);
SL_combo_rms_dB = NaN(length(files),1);
SL_combo_pp_dB = NaN(length(files),1);

for i = 1:length(files) %1:length(files)
    
    disp(i);    
    cd(files(i).name)
    
    files_mat = dir('*_RL_SL.mat');
    load(files_mat(1).name,'RL','SL','d');
    
    RL_rms_dB(i,:) = RL.rms_dB(3);
    RL_pp_dB(i,:) = RL.pp_dB(3);
    
    SL_model_rms_dB(i,:) = SL.model.rms_dB;
    SL_model_pp_dB(i,:) = SL.model.pp_dB;
    SL_combo_rms_dB(i,:) = SL.combo.rms_dB;
    SL_combo_pp_dB(i,:) = SL.combo.pp_dB;
    
    d_model(i,1) = d.model;
    d_combo(i,1) = d.combo;
    
    clearvars RL SL d
    
    cd ../
    
end

clear  i files* reindex range

save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\',...
    'typeUF1_combo_SCOT_564_to_937_RL_SL.mat'],...
    'RL*','SL*','d*');

clear
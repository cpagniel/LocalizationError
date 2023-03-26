%% An experimental approach for reducing error in passive acoustic localizations 
%% of individual marine animals using in situ sound playbacks in complex 
%% shallow-water environments
% Pagniello 2021

% This code combines .mat file outputed from
% Pagniello_2021_localization_TDOA_**_replicas_SCOT.m. Change file name and
% output manually as needed.

%% Load Log

Log = readtable('D:\Localization_Experiment\805572635\wav\Pagniello_etal_real_fish_log.xls'); % importing the log
Log = Log(:,1:8); % only care about these columns in the log

load('F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1\log_no564\log_no564_combo_SCOT_interp_5m_error_5m.mat','event');

start_sample = NaN(937-564,1);
%start_sample = NaN(length(event.start_time),1);

cnt = 1;
for i = 564:937 %1:length(event.start_time)
    start_sample(cnt) = floor(seconds(duration(event.start_time(i)-datetime(wavname2dnum_edit(Log{i,1}),'ConvertFrom','datenum'),'Format','s')).*48000); % event start sample
    cnt = cnt + 1;
end
clear i cnt

%% Combine Files within Folders

rootdir = 'F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1';
cd(rootdir);

files = dir('log*');
[~, reindex] = sort(str2double(regexp({files.name},'\d+','match','once')));
files = files(reindex);

range = 1:length(files);
files = files(range);
start_sample = start_sample(range);

xx = NaN(length(files),1);
yy = NaN(length(files),1);
zz = NaN(length(files),1);

MAX_C = cell(length(files),1);
MTDOA = cell(length(files),1);

DELTA_t = NaN(length(files),1);

fpeak = NaN(length(files),1);

BW3dB = NaN(length(files),1);
flow3dB = NaN(length(files),1);
fhigh3dB = NaN(length(files),1);

BW6dB = NaN(length(files),1);
flow6dB = NaN(length(files),1);
fhigh6dB = NaN(length(files),1);

BW10dB = NaN(length(files),1);
flow10dB = NaN(length(files),1);
fhigh10dB = NaN(length(files),1);

RL_rms_dB = NaN(length(files),1);
RL_pp_dB = NaN(length(files),1);

end_sample = NaN(length(files),1);
start_time = NaT(length(files),1);

for i = 1:length(files) %1:length(files)
    
    disp(i);    
    cd(files(i).name)
    
    files_mat = dir('*_combo_SCOT_interp_5m_error_5m.mat');
    load(files_mat(1).name,'event','position','mTDOA','max_C','PARAMS');
        
    xx(i) = position.X;
    yy(i) = position.Y;
    zz(i) = position.Z;
    
    MAX_C{i} = max_C.';
    MTDOA{i} = mTDOA.';
    
    DELTA_t(i) = PARAMS.dur;
    
    fpeak(i) = PARAMS.fpeak;
    
    BW3dB(i) = PARAMS.BW3dB;
    flow3dB(i) = PARAMS.flow3dB;
    fhigh3dB(i) = PARAMS.fhigh3dB;
    
    BW6dB(i) = PARAMS.BW6dB;
    flow6dB(i) = PARAMS.flow6dB;
    fhigh6dB(i) = PARAMS.fhigh6dB;
    
    BW10dB(i) = PARAMS.BW10dB;
    flow10dB(i) = PARAMS.flow10dB;
    fhigh10dB(i) = PARAMS.fhigh10dB;
    
    RL_rms_dB(i) = PARAMS.RL_rms_dB;
    RL_pp_dB(i) = PARAMS.RL_pp_dB;
        
    start_sample(i) = start_sample(i) + PARAMS.start_ind;
    end_sample(i) = start_sample(i) + PARAMS.end_ind;
    
    start_time(i) = seconds(start_sample(i)./48000)+datetime(wavname2dnum_edit(Log{i,1}),'ConvertFrom','datenum');
    
    clearvars position mTDOA max_C PARAMS
    
    cd ../
    
end

clear Log i files* reindex range

save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\',...
    'typeUF1_combo_SCOT_564_to_937_interp_5m_error_5m.mat'],...
    'xx','yy','zz',...
    'MAX_C','MTDOA',...
    'RL*','BW*','f*','DELTA_t',...
    'start*','end*');

clear
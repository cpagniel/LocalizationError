%% An experimental approach for reducing error in passive acoustic localizations 
%% of individual marine animals using in situ sound playbacks in complex 
%% shallow-water environments
% Pagniello 2021

% This code combines .mat file outputed from
% Pagniello_2021_localization_TDOA_**_replicas_SCOT.m. Change file name and
% output manually as needed.

%% Combine Files within Folders

rootdir = 'F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\0_interp_5';
cd(rootdir);

files = dir('event*');

event_no = NaN(length(files),1);

DELTA_x = NaN(length(files),1);
DELTA_y = NaN(length(files),1);
DELTA_z = NaN(length(files),1);
DELTA_xy = NaN(length(files),1);
DELTA_xyz = NaN(length(files),1);

x = NaN(length(files),1);
y = NaN(length(files),1);
z = NaN(length(files),1);

xx = NaN(length(files),1);
yy = NaN(length(files),1);
zz = NaN(length(files),1);

MAX_C = cell(length(files),1);
MTDOA = cell(length(files),1);

RL_rms_dB = cell(length(files),1);
RL_pp_dB = cell(length(files),1);

SNR_rms1 = cell(length(files),1);
SNR_rms_dB1 = cell(length(files),1);
SNR_rms2 = cell(length(files),1);
SNR_rms_dB2 = cell(length(files),1);

start_time = NaT(length(files),1);

c = NaN(length(files),1);

for i = 1:length(files)
    
    disp(i);    
    cd(files(i).name)
    
    files_mat = dir('*_data_SCOT.mat');
    load(files_mat(1).name,'event','error','position','mTDOA','max_C','PARAMS');
    
    event_no(i) = event.no;
    
    DELTA_x(i) = error.DELTA_x;
    DELTA_y(i) = error.DELTA_y;
    DELTA_z(i) = error.DELTA_z;
    DELTA_xy(i) = error.DELTA_xy;
    DELTA_xyz(i) = error.DELTA_xyz;
    
    x(i) = event.x;
    y(i) = event.y;
    z(i) = event.z;
    
    xx(i) = position.X;
    yy(i) = position.Y;
    zz(i) = position.Z;
    
    MAX_C{i} = max_C.';
    MTDOA{i} = mTDOA.';
    
    RL_rms_dB{i} = PARAMS.RL_rms_dB;
    RL_pp_dB{i} = PARAMS.RL_pp_dB;
    
    SNR_rms1{i} = PARAMS.SNR_rms1;
    SNR_rms_dB1{i} = PARAMS.SNR_rms_dB1;
    SNR_rms2{i} = PARAMS.SNR_rms2;
    SNR_rms_dB2{i} = PARAMS.SNR_rms_dB2;
    
    start_time(i) = event.start_time_corrected;
    
    c(i) = event.c;
    
    clearvars error event position mTDOA max_C PARAMS
    
    cd ../
    
end

save(['F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\',...
    'error_type0_data_SCOT_interp_5.mat'],...
    'DELTA*',...
    'x','y','z','xx','yy','zz','event_no',...
    'MAX_C','MTDOA',...
    'RL*','SNR*',...
    'start_time','c');

%clear
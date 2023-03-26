%% An experimental approach for reducing error in passive acoustic localizations 
%% of individual marine animals using in situ sound playbacks in complex 
%% shallow-water environments
% Pagniello 2021

% This code combines .mat file outputed from
% Pagniello_2021_localization_TDOA_**_replicas_SCOT.m. Change file name and
% output manually as needed.

%% Array Geometry

array.M = 4; % number of hydrophones
array.xr = zeros(array.M,1); array.yr = zeros(array.M,1); array.zr = zeros(array.M,1); % hydrophone positions in x, y, and z

array.GPS = load('F:\C_files\C_Spring_2021\Manuscript_Files\LocalizationError\code\Pagniello_2021_hydrophone_GPS_location.txt');

zone = utmzone(mean(array.GPS(:,1)),mean(array.GPS(:,2)));
utmstruct = defaultm('utm'); utmstruct.zone = zone; utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

% Note, center of array is position 1 for July 2018. Subtract array.xr(1)
% and array.yr(1) to center array at (0,0)
[array.x,array.y] = mfwdtran(utmstruct,array.GPS(:,1),array.GPS(:,2));
array.xr = array.x-array.x(1);
array.yr = array.y-array.y(1);
array.zr = distdim(array.GPS(:,3),'ft','m');

clear zone utmstruct

%% Combine Files within Folders

rootdir = 'F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq14';
cd(rootdir);

files = dir('log*');
[~, reindex] = sort(str2double(regexp({files.name},'\d+','match','once')));
files = files(reindex);

range = 1:length(files);
files = files(range);

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

d = NaN(length(files),1);

SL_rms_dB = NaN(length(files),1);
SL_pp_dB = NaN(length(files),1);

start_time = NaT(length(files),1);

grid_type = NaN(length(files),1);

for i = 1:length(files) %1:length(files)
    
    disp(i);    
    cd(files(i).name)
    
    files_mat = dir('*_combo_SCOT_interp_5m_error_5m.mat');
    load(files_mat(1).name,'data','event','position','mTDOA','max_C','PARAMS');
        
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
    
    d(i) = sqrt((xx(i)-array.xr(3)).^2 + (yy(i)-array.yr(3)).^2 + (zz(i)-array.zr(3)).^2);
    h = array.zr(3);
    
    if d(i) <= h
        SL_rms_dB(i) = RL_rms_dB(i) + 20*log10(h);
        SL_pp_dB(i) = RL_pp_dB(i) + 20*log10(h);
    else
        SL_rms_dB(i) = RL_rms_dB(i) + 20*log10(h) + 10*log10((d(i)-h)/h);
        SL_pp_dB(i) = RL_pp_dB(i) + 20*log10(h) + 10*log10((d(i)-h)/h);
    end

    start_time(i) = event.start_time(1108+i) + seconds(data.extend1) + seconds(PARAMS.start_ind./data.fs); %seconds(start_sample(i)./48000)+datetime(wavname2dnum_edit(Log{i,1}),'ConvertFrom','datenum');
    
    % For combo type only
    load('F:\C_files\C_Spring_2021\Real_Call_Localization\grid_type.mat');
    xgridno = find(-60:0.25:60 == position.X);
    ygridno = find(-60:0.25:60 == position.Y);
    if position.Z == 5
        zgridno = 1;
    elseif position.Z == 13
        zgridno = 2;
    end
    grid_type(i) = gridd.type(ygridno,xgridno,zgridno);
    
    clearvars data event position mTDOA max_C PARAMS *gridno
    
    cd ../
    
end

clear Log i files* reindex range array

save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\',...
    'typeUF1_combo_SCOT_1109_to_1527_interp_5m_error_5m.mat'],...
    'xx','yy','zz',...
    'MAX_C','MTDOA',...
    'RL*','SL*','BW*','f*','DELTA_t',...
    'start*','grid_type');

clear
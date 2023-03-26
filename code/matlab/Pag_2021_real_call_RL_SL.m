%% An experimental approach for quantifying and reducing error in passive acoustic localizations
%% of individual marine animals using in situ sound playbacks in complex
%% shallow-water environments
% Pagniello 2021

% Compute RL and SL between 245 and 1255 Hz.

%% Array Geometry

array.M = 4; % number of hydrophones
array.xr = zeros(array.M,1); array.yr = zeros(array.M,1); array.zr = zeros(array.M,1); % hydrophone positions in x, y, and z

array.GPS = load('F:\C_files\C_Spring_2021\Manuscript_Files\Localization\code\Pagniello_2021_hydrophone_GPS_location.txt');

zone = utmzone(mean(array.GPS(:,1)),mean(array.GPS(:,2)));
utmstruct = defaultm('utm'); utmstruct.zone = zone; utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

% Note, center of array is position 1 for July 2018. Subtract array.xr(1)
% and array.yr(1) to center array at (0,0)
[array.x,array.y] = mfwdtran(utmstruct,array.GPS(:,1),array.GPS(:,2));
array.xr = array.x-array.x(1);
array.yr = array.y-array.y(1);
array.zr = distdim(array.GPS(:,3),'ft','m');

% Hydrophone calibration
array.hsens = [-165.2 -165.5 -164.3 -165.5]; % hydrophone sensitivity for channels 1 to 4 (from calibration sheet)
array.tf = 20*log10(2^16/(1.5--1.5)); % transfer function for ST4300 DAQ

% Distance between each hydrophone
cnt = 1;
for m = 1:array.M-1
    for n = m+1:array.M
        array.dist(cnt) = sqrt((array.xr(m)-array.xr(n))^2 + ...
            (array.yr(m)-array.yr(n))^2 + (array.zr(m)-array.zr(n))^2);
        cnt = cnt + 1;
    end
end
clear m n cnt

array.maxTDOA = (array.dist./1516.0).'; % maximum TDOA based on array geometry and average sound speed

%% Loop Through All Logged Real Fish Calls

for gcnt = 1109:1527 %564:937 
    
    tic; % to calculate processing time for each event
    
    disp(gcnt);
    
    %% Import Data
    
    cd(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq14\log_no' num2str(gcnt)]);
%     cd(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq9\log_no' num2str(gcnt)]);
    
    load(['log_no' num2str(gcnt) '_combo_SCOT_interp_5m_error_5m.mat'],'PARAMS','data','position');
    DATA.dmean = data.dmean(PARAMS.start_ind:PARAMS.end_ind,:);
    pos.combo = position;
    FS = data.fs;
    
    clear data position PARAMS
    
    load(['log_no' num2str(gcnt) '_model_SCOT_interp_5m_error_5m.mat'],'position');
    pos.model = position;
    
    clear position
    
    %% Initial Signal Processsing
    
    bpFilt = designfilt('bandpassiir','FilterOrder',8,...
        'PassbandFrequency1',245,'PassbandFrequency2',1255, ...
        'PassbandRipple',1,'SampleRate',FS); % 8th-order Band Pass IIR
    for m = 1:array.M
        DATA.filt(:,m) = filter(bpFilt,DATA.dmean(:,m)); % filter data
    end
    
    clear bpFilt m
        
    %% Received Level
    
    RL.pp_dB = 20*log10(max(DATA.filt)-min(DATA.filt));
    RL.rms_dB = 20*log10(sqrt(mean(DATA.filt.^2)));
    
   %% Source Level
   
    h = array.zr(3); % for 1109:1527 & 564:937
       
    d.model = sqrt((pos.model.X-array.xr(3)).^2 + (pos.model.Y-array.yr(3)).^2 + (pos.model.Z-array.zr(3)).^2);   
    if d.model <= h
        SL.model.rms_dB = RL.rms_dB(3) + 20*log10(h);
        SL.model.pp_dB = RL.pp_dB(3) + 20*log10(h);
    else
        SL.model.rms_dB = RL.rms_dB(3) + 20*log10(h) + 10*log10((d.model-h)/h);
        SL.model.pp_dB = RL.pp_dB(3) + 20*log10(h) + 10*log10((d.model-h)/h);
    end
    
    d.combo = sqrt((pos.combo.X-array.xr(3)).^2 + (pos.combo.Y-array.yr(3)).^2 + (pos.combo.Z-array.zr(3)).^2);
    if d.combo <= h
        SL.combo.rms_dB = RL.rms_dB(3) + 20*log10(h);
        SL.combo.pp_dB = RL.pp_dB(3) + 20*log10(h);
    else
        SL.combo.rms_dB = RL.rms_dB(3) + 20*log10(h) + 10*log10((d.combo-h)/h);
        SL.combo.pp_dB = RL.pp_dB(3) + 20*log10(h) + 10*log10((d.combo-h)/h);
    end
    
    %% Save
    
    save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq14',...
        '\log_no' num2str(gcnt) '\log_no' num2str(gcnt) '_RL_SL.mat'],...
        'DATA','FS','pos','RL','SL','d','h');
    
%     save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1_seq9',...
%         '\log_no' num2str(gcnt) '\log_no' num2str(gcnt) '_RL_SL.mat'],...
%         'DATA','FS','pos','RL','SL','d','h');

    %% Clear
    
    clear DATA FS pos RL SL d h
    
end

clear

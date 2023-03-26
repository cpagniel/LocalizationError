%% An experimental approach for quantifying and reducing error in passive acoustic localizations
%% of individual marine animals using in situ sound playbacks in complex
%% shallow-water environments
% Pagniello 2021

% Grid-search algorithm combines modeled and data derived time-differences
% of arrival (TDOAs) from waveform smoothed coherence transform (SCOT) to
% localize real fish call.

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

%% Numerical Setup

% Search Grid
gridd.xmin = -60; gridd.ymin = -60; % grid limits (m)
gridd.xmax = 60; gridd.ymax = 60; % grid limits (m)
gridd.dx = 0.25; gridd.dy = 0.25; % resolution (m)
gridd.nx = 481; gridd.ny = 481; gridd.nz = 2; % number of grid points for x, y and z

[gridd.X,gridd.Y] = meshgrid(gridd.xmin:gridd.dx:gridd.xmax,gridd.ymin:gridd.dy:gridd.ymax);

%% Colormap for Ambiguity Surfaces

load('F:\C_files\C_Spring_2021\Localization_Playbacks\code\cmap_div.mat');

%% Load Call Log

Log = readtable('D:\Localization_Experiment\805572635\wav\Pagniello_etal_real_fish_log.xls'); % importing the log
Log = Log(:,1:8); % only care about these columns in the log

event.start_time = Log{:,5}; event.end_time = Log{:,6};
event.freq_max = Log{:,7}; event.freq_min = Log{:,8};

event.c = 1516.0; % harmonic average

%% Loop Through All Logged Real Fish Calls

for gcnt = 564:937 %1:length(event.start_time)
    
    tic; % to calculate processing time for each event
    
    disp(gcnt);
    
    %% Import .wav
    
    data.extend1 = 0; % import 0.01 seconds before logged call for SNR measurements
    data.extend2 = 0; % import 0.005 seconds after logged call for SNR measurements
    
    data.fs = getfield(audioinfo(char(Log{gcnt,1})),'SampleRate'); % sampling rate (Hz)
    
    data.start_sample = floor(seconds(duration(event.start_time(gcnt)-datetime(wavname2dnum_edit(Log{gcnt,1}),'ConvertFrom','datenum'),'Format','s')).*data.fs); % event start sample
    data.end_sample = ceil(seconds(duration(event.end_time(gcnt)-datetime(wavname2dnum_edit(Log{gcnt,1}),'ConvertFrom','datenum'),'Format','s')).*data.fs); % event end sample
    
    [data.raw,~] = audioread(char(Log{gcnt,1}),[data.start_sample-(data.extend1*data.fs) data.end_sample+(data.extend2*data.fs)],'native'); % event data raw audio from start to end of logged sound
    data.double = double(data.raw); % make data in double format
    data.daq = 10.^(-1.*(array.hsens+array.tf)/20).*data.double; % apply calibration
    
    data.dmean = bsxfun(@minus,data.daq,mean(data.daq)); % substract mean from data
    
    %% Initial Signal Processsing
    
    if event.freq_min(gcnt)-100 <= 0
        event.freq_min(gcnt) = 101;
    end
    
    bpFilt = designfilt('bandpassiir','FilterOrder',8,...
        'PassbandFrequency1',event.freq_min(gcnt)-100,'PassbandFrequency2',event.freq_max(gcnt)+100, ...
        'PassbandRipple',1,'SampleRate',data.fs); % 8th-order Band Pass IIR
    for m = 1:array.M
        data.filt(:,m) = filter(bpFilt,data.dmean(:,m)); % filter data
    end
    
    clear bpFilt m
    
    %% Create Output Directory
    
    warning('off');
    cd('F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1');
    mkdir(['log_no' num2str(gcnt)]);
    cd(['log_no' num2str(gcnt)]);
    
    %% Review Call and Calculate Call Parameters
    
    % Call Duration, Peak Frequency, Bandwdith and RL for All Channels
    
    Nfft = 512; Nfft2 = 8192; overlap = 0.9;
    
    for m = 1:array.M
        
        TEMP.start_ind(m) = find(cumsum(data.filt(:,m).^2) <= ...
            max(cumsum(data.filt(:,m).^2))*0.05,1,'last');
        TEMP.end_ind(m) = find(cumsum(data.filt(:,m).^2) <= ...
            max(cumsum(data.filt(:,m).^2))*0.95,1,'last');
        
        TEMP.dur(m) = (TEMP.end_ind(m) - TEMP.start_ind(m))./data.fs; % call duration
        
        
        [Pxx,f] = pwelch(data.filt(TEMP.start_ind(m):TEMP.end_ind(m),m),kaiser(Nfft,7.85),...
            round(Nfft*overlap),Nfft,data.fs);
        Pxx = Pxx(f <= 1200); f = f(f <= 1200); % remove high frequencies
        f_interp = 0:0.1:1200; Pxx_interp = interp1(f,Pxx,f_interp);
        
        [TEMP.BW3dB(m),TEMP.flow3dB(m),TEMP.fhigh3dB(m)] = powerbw(Pxx_interp,f_interp,[],3);
        [TEMP.BW6dB(m),TEMP.flow6dB(m),TEMP.fhigh6dB(m)] = powerbw(Pxx_interp,f_interp,[],6);
        [TEMP.BW10dB(m),TEMP.flow10dB(m),TEMP.fhigh10dB(m)] = powerbw(Pxx_interp,f_interp,[],10);
        
        if TEMP.flow10dB(m) == 0
            TEMP.flow10dB(m) = 1;
        end
        
        if TEMP.end_ind-TEMP.start_ind > Nfft2
            [Pxx,f] = pwelch(data.filt(TEMP.start_ind(m):TEMP.end_ind(m),m),kaiser(Nfft2,7.85),...
                round(Nfft2*overlap),Nfft2,data.fs);
            
            [~,ind] = max(Pxx);
            TEMP.fpeak(m) = f(ind);
        end
        
        bpFilt = designfilt('bandpassiir','FilterOrder',8,...
            'PassbandFrequency1',TEMP.flow10dB(m),'PassbandFrequency2',TEMP.fhigh10dB(m), ...
            'PassbandRipple',1,'SampleRate',data.fs); % 8th-order Band Pass IIR
        
        TEMP_data = filter(bpFilt,data.dmean(TEMP.start_ind(m):TEMP.end_ind(m),m)); % filter data
        
        TEMP.RL_pp_dB(m) = 20*log10(max(TEMP_data)-min(TEMP_data));
        TEMP.RL_rms_dB(m) = 20*log10(sqrt(mean(TEMP_data.^2)));
        
        clear bpFilt
        
    end
    
    clear Pxx* f* Nfft* overlap ind m TEMP_data
    
    % Save in PARAMS Parameters from Loudest Call
    
    ind_global = 2;
    % [~,ind] = max(TEMP.RL_rms_dB);
       
    PARAMS.ch = ind_global;
    if TEMP.start_ind(ind_global)/data.fs-0.02 >= 0
        PARAMS.start_ind = round(TEMP.start_ind(ind_global)-0.02.*data.fs);
    else
        PARAMS.start_ind = 1;
    end
    PARAMS.end_ind = round(TEMP.end_ind(ind_global)+0.02.*data.fs);
    PARAMS.dur = (PARAMS.end_ind-PARAMS.start_ind)./data.fs; %TEMP.dur(ind);
    
    PARAMS.fpeak = TEMP.fpeak(ind_global);
    
    PARAMS.BW3dB = TEMP.BW3dB(ind_global); PARAMS.flow3dB = TEMP.flow3dB(ind_global); PARAMS.fhigh3dB = TEMP.fhigh3dB(ind_global);
    PARAMS.BW6dB = TEMP.BW6dB(ind_global); PARAMS.flow6dB = TEMP.flow6dB(ind_global); PARAMS.fhigh6dB = TEMP.fhigh6dB(ind_global);
    PARAMS.BW10dB = TEMP.BW10dB(ind_global); PARAMS.flow10dB = TEMP.flow10dB(ind_global); PARAMS.fhigh10dB = TEMP.fhigh10dB(ind_global);
    
    PARAMS.RL_pp_dB = TEMP.RL_pp_dB(ind_global); PARAMS.RL_rms_dB = TEMP.RL_rms_dB(ind_global);
    
    if PARAMS.flow10dB == 0
        PARAMS.flow10dB = 1;
    end
    bpFilt = designfilt('bandpassiir','FilterOrder',8,...
        'PassbandFrequency1',PARAMS.flow10dB,'PassbandFrequency2',PARAMS.fhigh10dB, ...
        'PassbandRipple',1,'SampleRate',data.fs); % 8th-order Band Pass IIR
    
    for m = 1:array.M
        data.signal(:,m) = filter(bpFilt,data.dmean(PARAMS.start_ind:PARAMS.end_ind,m));
    end
    
    clear m bpFilt Nfft* overlap time
        
    %% Modelled Replica Fields
    
    if ~exist('cTDOA','var')
        
        array.maxTDOA = (array.dist./event.c).'; % maximum TDOA based on array geometry
        
        cTDOA = cell(gridd.ny,gridd.nx,gridd.nz);
        
        for i = 1:gridd.nx
            xg = gridd.xmin + (i-1)*gridd.dx;
            for j = 1:gridd.ny
                yg = gridd.ymin + (j-1)*gridd.dy;
                for k = 1:gridd.nz
                    if k == 1
                        zg = 5;
                    elseif k == 2
                        zg = 13;
                    end
                    
                    cnt = 1;
                    cTDOA_temp = zeros(nchoosek(array.M,2),1);
                    for m = 1:array.M-1
                        T_hat_m = norm([xg;yg;zg]-[array.xr(m,:);array.yr(m,:);array.zr(m,:)])/event.c;
                        for n = m+1:array.M
                            T_hat_n = norm([xg;yg;zg]-[array.xr(n,:);array.yr(n,:);array.zr(n,:)])/event.c;
                            
                            cTDOA_temp(cnt) = T_hat_n-T_hat_m;
                            cnt = cnt + 1;
                        end
                    end
                    
                    % remove maximum TDOA based on array geometry
                    ind = abs(cTDOA_temp) > array.maxTDOA;
                    cTDOA_temp(ind) = NaN;
                    
                    cTDOA{j,i,k} = cTDOA_temp;
                    clear cTDOA_temp ind
                    
                end
            end
        end
        gridd.xmax = xg; gridd.ymax = yg;
        
        clear i j k
        clear xg yg zg
        clear m n
        clear cnt
        clear T_hat_m T_hat_n
        clear ind
        
    end
    
    %% Measured Esimates of Time-Difference of Arrival using Waveform Smoothed Coherence Transform (SCOT)
    
    mTDOA = zeros(nchoosek(array.M,2),1);
    max_C = zeros(nchoosek(array.M,2),1);
    
    cnt = 1;
    for m = 1:array.M-1
        for n = m+1:array.M
            [C,lag] = SCOT(data.dmean(PARAMS.start_ind:PARAMS.end_ind,n),data.dmean(PARAMS.start_ind:PARAMS.end_ind,m),[PARAMS.flow10dB PARAMS.fhigh10dB],data.fs);
            C_bnd = C(lag./data.fs >= -array.maxTDOA(cnt) &  lag./data.fs <= array.maxTDOA(cnt));
            lag_bnd = lag(lag./data.fs >= -array.maxTDOA(cnt) &  lag./data.fs <= array.maxTDOA(cnt));
            [max_C(cnt),ind] = max(abs(C_bnd));
            lagDiff = lag_bnd(ind);
            mTDOA(cnt) = lagDiff/data.fs;
            
            cnt = cnt+1;
        end
    end
    clear cnt
    clear m n
    clear C lag
    clear ind
    clear lagDiff
    
    %% Maximum Likelihood Estimate of the Source Position
    
    LS = zeros(gridd.ny,gridd.nx,gridd.nz);
    
    for xx = 1:gridd.nx
        for yy = 1:gridd.ny
            for zz = 1:gridd.nz
                LS(yy,xx,zz) = nansum((mTDOA-cTDOA{yy,xx,zz}).^2);
            end
        end
    end
    clear xx yy zz
    
    ind = LS == 0; LS(ind) = NaN;
    clear ind
    
    [~,ind] = min(LS(:));
    [ygridno,xgridno,zgridno] = ind2sub(size(LS),ind);
    clear ind
    
    position.X = (xgridno-1)*gridd.dx + gridd.xmin;
    position.Y = (ygridno-1)*gridd.dy + gridd.ymin;
    if zgridno == 1
        position.Z = 5;
    elseif zgridno == 2
        position.Z = 13;
    end
    
    elapsed_time = toc;
    
    %% Save
    
    save(['F:\C_files\C_Spring_2021\Real_Call_Localization\Localization_Output\UF1',...
        '\log_no' num2str(gcnt) '\log_no' num2str(gcnt) '_model_SCOT_interp_5m_error_5m.mat'],...
        'cTDOA','mTDOA','max_C','LS','position','event','data','elapsed_time','TEMP','PARAMS');
    
    %% Plot LS
    
    LS_norm = (LS - min(LS(:)))/(max(LS(:))-min(LS(:)));
    
    fig = figure('Visible','off');
    fig.Position = [279 49 899 768];
    
    pcolor(gridd.xmin:gridd.dx:gridd.xmax,gridd.ymin:gridd.dy:gridd.ymax,LS_norm(:,:,zgridno));
    shading flat; colormap(flipud(cmap_div)); colorbar; caxis([0 1]);
    
    hold on
    
    plot([array.xr(1) array.xr(2)],[array.yr(1) array.yr(2)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    plot([array.xr(1) array.xr(3)],[array.yr(1) array.yr(3)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    plot([array.xr(1) array.xr(4)],[array.yr(1) array.yr(4)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    
    hold on
    
    plot(position.X,position.Y,'ro','MarkerFaceColor','r','MarkerSize',4);
    
    hold on
    
    title(['log no ' num2str(gcnt)])
    xlabel('Eastings (m)'); ylabel('Northings (m)');
    axis square
    xlim([gridd.xmin gridd.xmax]); ylim([gridd.ymin gridd.ymax])
    
    saveas(gcf,['log_no' num2str(gcnt) '_LS_model_SCOT_interp_5m_error_5m.fig']);
    saveas(gcf,['log_no' num2str(gcnt) '_LS_model_SCOT_interp_5m_error_5m.png']);
    
    close all
    
    clear fig LS_norm
    
    
    %% Clear
    
    clear mTDOA
    clear max_C
    clear LS
    clear position
    clear elapsed_time
    clear xgridno ygridno zgridno
    clear data
    
end

clear

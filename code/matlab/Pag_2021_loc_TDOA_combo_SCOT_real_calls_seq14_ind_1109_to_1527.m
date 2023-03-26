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

%% Load Playback Experimental Log

[LOG_num,LOG_txt] = xlsread('F:\C_files\C_Spring_2021\Localization_Playbacks\code\Pagniello_2021_experiment_data_log.xlsx','Sheet1','A2:K8163','basic');

% Total Number of Playbacks
event = 1:size(LOG_num,1); event = event.';
GPS_lat = LOG_num(:,10); GPS_lon = LOG_num(:,11); z = LOG_num(:,8);

% Only Keep Playbacks of Fish Calls
ind = LOG_num(:,7) == 0;
event = event(ind);
GPS_lat = GPS_lat(ind); GPS_lon = GPS_lon(ind); z = z(ind);
clear ind

% Only Keep Playbacks at 5 +/- 1 m or 13 - 1 m
round_z = round(z);
round_z(round_z == 4) = 5; round_z(round_z == 6) = 5;
ind5 = round_z == 5; % "shallow" playbacks
round_z(round_z == 12) = 13;
ind13 = round_z == 13; % "deep" playbacks
ind = (ind5 + ind13) == 1; % "shallow" & "deep" playbacks

% Meta Data of Playbacks of Fish Calls
event_no = event(ind); % event number
GPS_lat = GPS_lat(ind); GPS_lon = GPS_lon(ind); z_play = z(ind); % GPS & Ambient-Pressure Position

% Convert GPS Position to X & Y
[x,y] = mfwdtran(utmstruct,GPS_lat,GPS_lon);
x_play = x-array.x(1);
y_play = y-array.y(1);

clear LOG* event ind* GPS* x y z round_z

%% Load Call Log

Log = readtable('D:\Localization_Experiment\805572635\wav\Pagniello_etal_real_fish_log.xls'); % importing the log
Log = Log(:,1:8); % only care about these columns in the log

event.start_time = Log{:,5}; event.end_time = Log{:,6};
event.freq_max = Log{:,7}; event.freq_min = Log{:,8};

event.c = 1516.0; % harmonic average

%% Loop Through All Logged Real Fish Calls

for gcnt = 1109:1527 %1:length(event.start_time)
    
    tic; % to calculate processing time for each event
    
    disp(gcnt);
    
    %% Import .wav
    
    data.extend1 = 0.05; % import 0.01 seconds before logged call for SNR measurements
    data.extend2 = -0.05; % import 0.005 seconds after logged call for SNR measurements
    
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
    time = (0:length(data.filt(:,1))-1)./data.fs;
    set(0,'DefaultFigureWindowStyle','normal')
    set(0,'DefaultFigureVisible','on')
    
    f1 = figure(1); set(f1,'Position',[1619 -29 650 650],'Visible','off');
    f2 = figure(2); set(f2,'Position',[2273 -29 650 650],'Visible','off');
    f3 = figure(3); set(f3,'Position',[2883 -29 650 650],'Visible','off');
    
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
        
        figure(1); subplot(2,2,m);
        [~,F,T,P] = spectrogram(data.filt(:,m),kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,data.fs,'yaxis');
        imagesc(T,F,10*log10(P))
        set(gca,'YDir','normal'); ylim([0 1500]);
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        colormap(jet); h = colorbar; caxis([60 90]); ylabel(h,{'Spectral Density (dB re 1 \muPa^2/Hz)'});
        axis square
        
        figure(2); subplot(2,2,m);
        plot(time,data.filt(:,m)/10^6,'k-')
        title({['Call Duration = ' num2str(TEMP.dur(m)) ' s'];['RL = ' num2str(TEMP.RL_pp_dB(m)) ' dB re 1 \muPa pp']});
        xlabel('Time (s)'); ylabel('Amplitude (Pa)');
        axis square
        
        figure(3); subplot(2,2,m);
        plot(f,10*log10(Pxx),'k-');
        xlabel('Frequency (Hz)'); ylabel({'Power Spectral Density';'(dB re 1 \muPa^2/Hz)'});
        xlim([0 1500]);
        
    end
    
    clear Pxx* f* Nfft* overlap ind m TEMP_data
    
    % Save in PARAMS Parameters from Channel 2 (Loudest and Cleaneast call)
    
    ind = 2;
    
    for m = 1:array.M
        
        figure(1); subplot(2,2,m);
        hold on
        plot([TEMP.start_ind(ind) TEMP.start_ind(ind)]./data.fs,[min(ylim) max(ylim)],'w--','LineWidth',2);
        plot([TEMP.end_ind(ind) TEMP.end_ind(ind)]./data.fs,[min(ylim) max(ylim)],'w--','LineWidth',2);
        plot([min(xlim) max(xlim)],[TEMP.flow10dB(ind) TEMP.flow10dB(ind)],'w--','LineWidth',2);
        plot([min(xlim) max(xlim)],[TEMP.fhigh10dB(ind) TEMP.fhigh10dB(ind)],'w--','LineWidth',2);
        
        figure(2); subplot(2,2,m);
        hold on
        plot(time(TEMP.start_ind(ind):TEMP.end_ind(ind)),data.filt(TEMP.start_ind(ind):TEMP.end_ind(ind),m)/10^6,'r-')
        
        figure(3); subplot(2,2,m);
        hold on
        plot([TEMP.fpeak(ind) TEMP.fpeak(ind)],[min(ylim) max(ylim)],'k--','LineWidth',2);
        plot([TEMP.flow10dB(ind) TEMP.flow10dB(ind)],[min(ylim) max(ylim)],'k-','LineWidth',2);
        plot([TEMP.fhigh10dB(ind) TEMP.fhigh10dB(ind)],[min(ylim) max(ylim)],'k-','LineWidth',2);
        
    end
    
    saveas(figure(1),['log_no' num2str(gcnt) '_spec.fig']); saveas(figure(1),['log_no' num2str(gcnt) '_spec.png']);
    saveas(figure(2),['log_no' num2str(gcnt) '_time.fig']); saveas(figure(2),['log_no' num2str(gcnt) '_time.png']);
    saveas(figure(3),['log_no' num2str(gcnt) '_freq.fig']); saveas(figure(3),['log_no' num2str(gcnt) '_freq.png']);
    
    clear f1 f2 f3
    
    PARAMS.ch = ind;
    PARAMS.start_ind = TEMP.start_ind(ind);
    PARAMS.end_ind = TEMP.end_ind(ind);
    PARAMS.dur = TEMP.dur(ind);
    
    PARAMS.fpeak = TEMP.fpeak(ind);
    
    PARAMS.BW3dB = TEMP.BW3dB(ind); PARAMS.flow3dB = TEMP.flow3dB(ind); PARAMS.fhigh3dB = TEMP.fhigh3dB(ind);
    PARAMS.BW6dB = TEMP.BW6dB(ind); PARAMS.flow6dB = TEMP.flow6dB(ind); PARAMS.fhigh6dB = TEMP.fhigh6dB(ind);
    PARAMS.BW10dB = TEMP.BW10dB(ind); PARAMS.flow10dB = TEMP.flow10dB(ind); PARAMS.fhigh10dB = TEMP.fhigh10dB(ind);
    
    PARAMS.RL_pp_dB = TEMP.RL_pp_dB(ind); PARAMS.RL_rms_dB = TEMP.RL_rms_dB(ind);
    
    if PARAMS.flow10dB == 0
        PARAMS.flow10dB = 1;
    end
    bpFilt = designfilt('bandpassiir','FilterOrder',8,...
        'PassbandFrequency1',PARAMS.flow10dB,'PassbandFrequency2',PARAMS.fhigh10dB, ...
        'PassbandRipple',1,'SampleRate',data.fs); % 8th-order Band Pass IIR
    
    for m = 1:array.M
        data.signal(:,m) = filter(bpFilt,data.dmean(PARAMS.start_ind:PARAMS.end_ind,m));
    end
    
    clear ind m bpFilt Nfft* overlap time
    
    %% Modelled Replica Fields
    
    if ~exist('cTDOA','var')
        
        array.maxTDOA = (array.dist./event.c).'; % maximum TDOA based on array geometry
        
        cTDOA_model = cell(gridd.ny,gridd.nx,gridd.nz);
        
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
                    
                    cTDOA_model{j,i,k} = cTDOA_temp;
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
        clear ind*
        
        %% Data Derived Replica Fields
        
        % Create cTDOA Vectors At Points Within the Search Grid
        cTDOA_data = cell(gridd.ny,gridd.nx,gridd.nz);
        
        % Compute Surface for Each Depth
        for zz = 1:gridd.nz
            
            % See Pagniello_2021_TDOAs_SCOT_for_data_replicas.m for TDOA
            % calculations.
            load('F:\C_files\C_Spring_2021\Localization_Playbacks\mTDOA_SCOT_for_type0.mat','mTDOA','x','y','z');
            
            % To eliminate data replicas whose DELTA_xy are greater than 5m
            load('F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\error_type0_data_SCOT.mat','start_time','DELTA_xy');
            
            [~,ind] = sort(start_time);
            DELTA_xy = DELTA_xy(ind);
            clear ind start_time
            
            x(DELTA_xy > 5) = []; y(DELTA_xy > 5) = []; z(DELTA_xy > 5) = [];
            mTDOA(DELTA_xy > 5,:) = [];
            
            % Separate mTDOA to make two depth surfaces
            if zz == 1 % z = 5 m, shallow
                ind = z <= 8;
                x = x(ind); y = y(ind);
                mTDOA = mTDOA(ind,:);
            elseif zz == 2 % z = 13 m, deep
                ind = z > 8;
                x = x(ind); y = y(ind);
                mTDOA = mTDOA(ind,:);
            end
            clear ind
            
            % Use scatteredInterpolant to Create Replica Vector Surface From Playbacks
            % natural interpolation, no extrapolation
            F_12 = scatteredInterpolant(x,y,mTDOA(:,1),'natural','none');
            F_13 = scatteredInterpolant(x,y,mTDOA(:,2),'natural','none');
            F_14 = scatteredInterpolant(x,y,mTDOA(:,3),'natural','none');
            F_23 = scatteredInterpolant(x,y,mTDOA(:,4),'natural','none');
            F_24 = scatteredInterpolant(x,y,mTDOA(:,5),'natural','none');
            F_34 = scatteredInterpolant(x,y,mTDOA(:,6),'natural','none');
            
            % Create cTDOA Vectors At Points Within the Search Grid
            for i = 1:gridd.nx
                for j = 1:gridd.ny
                    cTDOA_temp = [F_12(gridd.X(j,i),gridd.Y(j,i)); F_13(gridd.X(j,i),gridd.Y(j,i));...
                        F_14(gridd.X(j,i),gridd.Y(j,i)); F_23(gridd.X(j,i),gridd.Y(j,i));...
                        F_24(gridd.X(j,i),gridd.Y(j,i)); F_34(gridd.X(j,i),gridd.Y(j,i))];
                    
                    % remove maximum TDOA based on array geometry
                    ind = abs(cTDOA_temp) > array.maxTDOA;
                    cTDOA_temp(ind) = NaN;
                    
                    cTDOA_data{j,i,zz} = cTDOA_temp;
                    clear cTDOA_temp ind
                end
            end
            clear i j
            clear temp
        end
        
        clear zz
        clear x y
        clear mTDOA
        clear F*
        clear DELTA_xy
        clear ind
        
        %% Combined Replica Field
        
        cTDOA = cTDOA_model;
        gridd.type = zeros([size(gridd.X) 2]);
        
        % Is Grid Point Within 5 m of a Playback?
        for i = 1:size(gridd.X,2)
            for j = 1:size(gridd.X,1)
                for k = 1:2
                    if k == 1 % z = 5 m, shallow
                        ind = z <= 8;
                        xi = x_play(ind); yi = y_play(ind);
                    elseif k == 2 % z = 13 m, deep
                        ind = z > 8;
                        xi = x_play(ind); yi = y_play(ind);
                    end
                    
                    ind = sum(sqrt((gridd.X(1,i)-xi).^2 + (gridd.Y(j,1)-yi).^2) <= 5);
                    if ind >= 1
                        if sum(~isnan(cTDOA_data{j,i,k})) ~= 0
                            gridd.type(j,i,k) = 1;
                            cTDOA{j,i,k} = cTDOA_data{j,i,k};
                        end
                    end
                end
            end
        end
        clear xi yi i j k ind
    end
    
    clear cTDOA_data cTDOA_model z *play
    
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
        '\log_no' num2str(gcnt) '\log_no' num2str(gcnt) '_combo_SCOT_interp_5m_error_5m.mat'],...
        'cTDOA','mTDOA','max_C','LS','position','event','data','elapsed_time','TEMP','PARAMS');
    
    %% Plot LS
    
    LS_norm = (LS - min(LS(:)))/(max(LS(:))-min(LS(:)));
    
    fig = figure('Visible','off');
    fig.Position = [279 49 899 768];
    
    pcolor(gridd.xmin:gridd.dx:gridd.xmax,gridd.ymin:gridd.dy:gridd.ymax,LS_norm(:,:,zgridno));
    shading flat; colormap(flipud(cmap_div)); h = colorbar; caxis([0 1]);
    
    hold on
    
    plot([array.xr(1) array.xr(2)],[array.yr(1) array.yr(2)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    plot([array.xr(1) array.xr(3)],[array.yr(1) array.yr(3)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    plot([array.xr(1) array.xr(4)],[array.yr(1) array.yr(4)],'-dk','LineWidth',1,'MarkerFaceColor','k');
    
    hold on
    
    if gridd.type(ygridno,xgridno,zgridno) == 1
        p1 = plot(position.X,position.Y,'go','MarkerFaceColor','g','MarkerSize',4);
    else
        p1 = plot(position.X,position.Y,'ro','MarkerFaceColor','r','MarkerSize',4);
    end
    
    hold on
    
    title(['log no ' num2str(gcnt)])
    xlabel('Eastings (m)'); ylabel('Northings (m)');
    axis square
    xlim([gridd.xmin gridd.xmax]); ylim([gridd.ymin gridd.ymax])
    
    saveas(gcf,['log_no' num2str(gcnt) '_LS_combo_SCOT_interp_5m_error_5m.fig']);
    saveas(gcf,['log_no' num2str(gcnt) '_LS_combo_SCOT_interp_5m_error_5m.png']);
    
    close all
    
    clear fig h p1 LS_norm
    
    
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

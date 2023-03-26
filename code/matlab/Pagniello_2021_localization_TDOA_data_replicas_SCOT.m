%% An experimental approach for reducing error in passive acoustic localizations 
%% of individual marine animals using in situ sound playbacks in complex 
%% shallow-water environments
% Pagniello 2021

% Grid-search algorithm uses data derived time-differences 
% of arrival (TDOAs) from waveform smoothed coherence transform (SCOT) to 
% localize fish call playbacks.

%% Array Geometry

array.M = 4; % number of hydrophones
array.xr = zeros(array.M,1); array.yr = zeros(array.M,1); array.zr = zeros(array.M,1); % hydrophone positions in x, y, and z

array.GPS = load('F:\C_files\C_Spring_2021\Localization_Playbacks\code\Pagniello_2021_hydrophone_GPS_location.txt');

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

%% Load Experimental Log

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
GPS_lat = GPS_lat(ind); GPS_lon = GPS_lon(ind); z = z(ind); % GPS & Ambient-Pressure Position

% Convert GPS Position to X & Y
[x,y] = mfwdtran(utmstruct,GPS_lat,GPS_lon);
xx = x-array.x(1);
yy = y-array.y(1);

clear LOG* event ind* GPS* x y

%% Loop Through All Events 

for gcnt = 1:length(event_no)
    
    tic; % to calculate processing time for each event
    
    disp(event_no(gcnt));
    
    %% Select Experimental Event
    
    event.no = event_no(gcnt); % event number
    
    [event_num,event_txt,~] = xlsread('F:\C_files\C_Spring_2021\Localization_Playbacks\code\Pagniello_2021_experiment_data_log.xlsx','Sheet1',['A' num2str(event_no(gcnt)+1) ':K' num2str(event_no(gcnt)+1)],'basic');
    
    event.file = event_txt{1}; % file location and name
    event.start_time = datetime(event_num(3),'ConvertFrom','excel'); % event start time
    event.start_time_corrected = datetime(event_num(4),'ConvertFrom','excel'); % corrected event start time
    event.end_time = datetime(event_num(5),'ConvertFrom','excel'); % event end time
    event.end_time_corrected = datetime(event_num(6),'ConvertFrom','excel'); % corrected event end time
    
    event.freqID = event_num(7); % frequency (Hz)
    if event.freqID == 0 % fish sound
        event.freq = 250:1:1250;
    elseif event.freqID == 999 % LFM sweep
        event.freq =  200:1:2000;
    elseif event.freqID == 111 % combination of tones
        event.freq = [400, 450, 515, 570, 630, 685];
    else
        event.freq = event.freqID;
    end
    
    event.z = event_num(8)+0.29; % depth (m)
    event.c = 1516; % harmonic average of speed of sound (m/s) event_num(9);
    event.GPS_lat = event_num(10); % GPS latitude (deg N)
    event.GPS_lon = event_num(11); % GPS longitude (deg W)
    
    [event.x,event.y] = mfwdtran(utmstruct,event.GPS_lat,event.GPS_lon);
    event.x = event.x-array.x(1);
    event.y = event.y-array.y(1);
    
    clear event_num event_txt
    
    %% Import .wav
    
    data.extend = 0.175; % import 0.175 seconds before and after logged call for SNR measurements
    
    data.fs = getfield(audioinfo(char(event.file)),'SampleRate'); % sampling rate (Hz)
    
    data.start_sample = floor(seconds(duration(event.start_time-datetime(wavname2dnum_edit(event.file),'ConvertFrom','datenum'),'Format','s')).*data.fs); % event start sample
    data.end_sample = ceil(seconds(duration(event.end_time-datetime(wavname2dnum_edit(event.file),'ConvertFrom','datenum'),'Format','s')).*data.fs); % event end sample
    
    [data.raw,~] = audioread(char(event.file),[data.start_sample-(data.extend*data.fs) data.end_sample+(data.extend*data.fs)],'native'); % event data raw audio from start to end of logged sound
    data.double = double(data.raw); % make data in double format
    data.daq = 10.^(-1.*(array.hsens+array.tf)/20).*data.double; % apply calibration
    
    data.dmean = bsxfun(@minus,data.daq,mean(data.daq)); % substract mean from data
    
    %% Signal Processsing
    
    bpFilt = designfilt('bandpassiir','FilterOrder',8,...
        'PassbandFrequency1',event.freq(1)-5,'PassbandFrequency2',event.freq(end)+5, ...
        'PassbandRipple',1,'SampleRate',data.fs); % 8th-order Band Pass IIR
    for m = 1:array.M
        data.filt(:,m) = filter(bpFilt,data.dmean(:,m)); % filter data
    end
    
    clear bpFilt m
    
    %% Create Output Directory
    
    warning('off');
    cd(['F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\' num2str(event.freqID)]);
    mkdir(['event_no' num2str(event_no(gcnt))]);
    cd(['event_no' num2str(event_no(gcnt))]);
    
    %% Calculate SNR
    
    data.signal = data.filt((data.extend*data.fs):end-(data.extend*data.fs),:); % fish call playback
    
    data.noise1 = data.filt(1:(data.extend*data.fs),:); % noise before fish call
    data.noise2 = data.filt(end-(data.extend*data.fs):end,:); % noise after fish call
    
    for m = 1:array.M
        
        % Noise Level Before Fish Call (root mean square)
        PARAMS.NL_rms1(m) = sqrt(mean(data.noise1(:,m).^2));
        PARAMS.NL_rms_dB1(m) = 20*log10(PARAMS.NL_rms1(m));
        
        % Noise Level After Fish Call (root mean square)
        PARAMS.NL_rms2(m) = sqrt(mean(data.noise2(:,m).^2));
        PARAMS.NL_rms_dB2(m) = 20*log10(PARAMS.NL_rms2(m));
        
        % Peak to Peak & RMS Received Level of Fish Call
        PARAMS.RL_pp_dB(m) = 20*log10(max(data.signal(:,m))-min(data.signal(:,m)));
        PARAMS.RL_rms(m) = sqrt(mean(data.signal(:,m).^2));
        PARAMS.RL_rms_dB(m) = 20*log10(PARAMS.RL_rms(m));
        
        % SNR of Fish Call to Noise Before Fish Call
        PARAMS.SNR_rms1(m) = PARAMS.RL_rms(m)/PARAMS.NL_rms1(m);
        PARAMS.SNR_rms_dB1(m) = 20*log10(PARAMS.SNR_rms1(m));
        
        % SNR of Fish Call to Noise After Fish Call
        PARAMS.SNR_rms2(m) = PARAMS.RL_rms(m)/PARAMS.NL_rms2(m);
        PARAMS.SNR_rms_dB2(m) = 20*log10(PARAMS.SNR_rms2(m));
        
    end
    clear m
    clear ind test
    
    %% Data Derived Replica Fields
    
    % Maximum TDOA Based on Array Geometry
    array.maxTDOA = (array.dist./event.c).'; 
    
    % Create cTDOA Vectors At Points Within the Search Grid
    cTDOA = cell(gridd.ny,gridd.nx,gridd.nz);
        
    % Compute Surface for Each Depth
    for zz = 1:gridd.nz
        
        % See Pagniello_2021_TDOAs_SCOT_for_data_replicas.m for TDOA
        % calculations.
        load('F:\C_files\C_Spring_2021\Localization_Playbacks\mTDOA_SCOT_for_type0.mat','mTDOA','x','y','z');
        
        % Remove Current Event
        mTDOA(gcnt,:) = []; x(gcnt) = []; y(gcnt) = []; z(gcnt) = [];
        
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
                    
                cTDOA{j,i,zz} = cTDOA_temp;
                clear cTDOA_temp ind
            end
        end
        clear i j
    end
    
    clear xx yy zz
    clear x y 
    clear mTDOA
    
     %% Measured Esimates of Time-Difference of Arrival using Waveform Smoothed Coherence Transform (SCOT)
    
    mTDOA = zeros(nchoosek(array.M,2),1);
    max_C = zeros(nchoosek(array.M,2),1);
    
    cnt = 1;
    for m = 1:array.M-1
        for n = m+1:array.M
            [C,lag] = SCOT(data.dmean((data.extend*data.fs):end-(data.extend*data.fs),n),...
                data.dmean((data.extend*data.fs):end-(data.extend*data.fs),m),event.freq,data.fs);
            % [C,lag] = SCOT(data.signal(:,n),data.signal(:,m),event.freq,data.fs);
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
    clear C* lag*
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
    
    %% Offset from GPS (i.e., error)
    
    error.DELTA_x = position.X - event.x;
    error.DELTA_y = position.Y - event.y;
    error.DELTA_z = position.Z - event.z;
    
    error.DELTA_xy = sqrt((position.X-event.x)^2 + (position.Y-event.y)^2);
    error.DELTA_xyz = sqrt((position.X-event.x)^2 + (position.Y-event.y)^2 + (position.Z-event.z)^2);
    
    elapsed_time = toc;
    
    %% Save
    
    save(['F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\',...
        num2str(event.freqID) '\event_no' num2str(event_no(gcnt)) '\event_no' num2str(event_no(gcnt)) '_data_SCOT.mat'],...
        'cTDOA','F*','mTDOA','max_C','LS','position','error','event','data','elapsed_time','PARAMS');
    
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
    
    p1 = plot(position.X,position.Y,'ro','MarkerFaceColor','r','MarkerSize',4);
    
    hold on
    
    p2 = plot(event.x,event.y,'go','MarkerFaceColor','g','MarkerSize',4);
    
    title(['event no ' num2str(event_no(gcnt))])
    xlabel('Eastings (m)'); ylabel('Northings (m)');
    axis square
    xlim([gridd.xmin gridd.xmax]); ylim([gridd.ymin gridd.ymax])
    
    saveas(gcf,['event_no' num2str(event_no(gcnt)) '_LS_data_SCOT.fig']);
    saveas(gcf,['event_no' num2str(event_no(gcnt)) '_LS_data_SCOT.png']);
    
    close all
    
    %% Clear
    
    clear cTDOA
    clear mTDOA
    clear max_C
    clear LS
    clear position
    clear error
    clear elapsed_time
    clear xgridno ygridno zgridno
    clear data
    clear event
 
end

clear

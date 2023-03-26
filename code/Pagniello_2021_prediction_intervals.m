%% Passive acoustic localization of individual marine animals using
%% in situ sound playbacks in complex shallow-water environments
% Pagniello 2021

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

%% Figure

figure;
sub_cnt = 1;

for i = 1:2
    
    %% Load Localization Data
    
    if i == 1
        load('F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\error_type0_model_SCOT.mat','x','y','z','xx','yy','zz');
        col = [0 0 0]; % black
    elseif i == 2
        load('F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\error_type0_data_SCOT.mat','x','y','z','xx','yy','zz');
        col = [1 0 0]; % red
    end
    
    
    %% Plot 
       
    subplot(2,3,sub_cnt);
    
    % Prediction Intervals Based on 95% Percentiles
    
    edges = -60:0.25:60;
    ubins = unique(xx);
    int = zeros(length(edges),2);
    
    for ii = 1:length(edges)
        for jj = 1:length(ubins)
            if edges(ii) == ubins(jj)
                int(ii,:) = [prctile(x(xx == ubins(jj)),2.5) prctile(x(xx == ubins(jj)),97.5)];
            end
        end
    end
    
    fill([edges.'; 0; flipud(edges.'); 0],[int(:,1); 0; flipud(int(:,2)); 0],col,'EdgeColor','none','FaceAlpha',0.5)
    
    hold on
    
    % XY Data
    
    plot(-60:60,-60:60,'k--','LineWidth',2);
    hold on
    plot(xx,x,'.','Color',col,'MarkerSize',3);
    
    % Axes 
    
    set(gca,'FontSize',18);
    xlim([-60 60]); ylim([-60 60]); 
    xticks([-60 -30 0 30 60]); yticks([-60 -30 0 30 60]);
    axis square
    grid on
    grid minor
    if i == 1
        xticklabels([]);
    end
    if i == 2
        xlabel({'Localization-Derived'; 'Estimated X Position (m)'},'FontSize',24); 
        lab = ylabel({'GPS-Derived'; 'True X Position (m)'},'FontSize',24);
        
        p = get(lab,'Position'); p(1) = -80; p(2) = 65;
        set(lab,'Position',p);
        
        a = get(gca,'Position'); a(2) = 0.19;
        set(gca,'Position',a)
    end
    set(gca,'LineWidth',1)
       
    subplot(2,3,sub_cnt+1);
    
    % Prediction Intervals Based on 95% Percentiles
    
    edges = -60:0.25:60;
    ubins = unique(yy);
    int = zeros(length(edges),2);
    
    for ii = 1:length(edges)
        for jj = 1:length(ubins)
            if edges(ii) == ubins(jj)
                int(ii,:) = [prctile(y(yy == ubins(jj)),2.5) prctile(y(yy == ubins(jj)),97.5)];
            end
        end
    end
    
    fill([edges.'; 0; flipud(edges.'); 0],[int(:,1); 0; flipud(int(:,2)); 0],col,'EdgeColor','none','FaceAlpha',0.5)
    
    hold on
    
    % XY Data
    
    plot(-60:60,-60:60,'k--','LineWidth',2);
    hold on
    plot(yy,y,'.','Color',col,'MarkerSize',3);

    % Axes
    
    set(gca,'FontSize',18);
    xlim([-60 60]); ylim([-60 60]);
    xticks([-60 -30 0 30 60]); yticks([-60 -30 0 30 60]);
    axis square
    grid on
    grid minor
    if i == 1
        xticklabels([]);
    end
    if i == 2
        xlabel({'Localization-Derived'; 'Estimated Y Position (m)'},'FontSize',24); 
        lab = ylabel({'GPS-Derived'; 'True Y Position (m)'},'FontSize',24);
        
        p = get(lab,'Position'); p(1) = -80; p(2) = 65;
        set(lab,'Position',p);
        
        a = get(gca,'Position'); a(2) = 0.19;
        set(gca,'Position',a)
    end
    set(gca,'LineWidth',1)
    
    subplot(2,3,sub_cnt+2);
    
    % Prediction Intervals Based on 95% Percentiles
    
    edges = -15:0.25:15;
    ubins = unique(zz);
    int = zeros(length(edges),2);
    
    for ii = 1:length(edges)
        for jj = 1:length(ubins)
            if edges(ii) == ubins(jj)
                int(ii,:) = [prctile(z(zz == ubins(jj)),2.5) prctile(z(zz == ubins(jj)),97.5)];
            end
        end
    end
    
    fill([edges.'; 0; flipud(edges.'); 0],[int(:,1); 0; flipud(int(:,2)); 0],col,'EdgeColor','none','FaceAlpha',0.5)
    
    hold on
    
    % XY Data
    
    plot(0:15,0:15,'k--','LineWidth',2);
    hold on
    plot(zz,z,'.','Color',col,'MarkerSize',3);
    
    % Axes
    
    set(gca,'FontSize',18);
    xlim([0 15]); ylim([0 15]);
    xticks([0 5 10 15]); yticks([0 5 10 15]);
    axis square
    grid on
    grid minor
    if i == 1
        xticklabels([]);
    end
    if i == 2
        xlabel({'Localization-Derived'; 'Estimated Z Position (m)'},'FontSize',24); 
        lab = ylabel({'Ambient Pressure-Derived'; 'True Z Position (m)'},'FontSize',24);
        
        p = get(lab,'Position'); p(1) = -2; p(2) = 16;
        set(lab,'Position',p);
        
        a = get(gca,'Position'); a(2) = 0.19;
        set(gca,'Position',a)
    end
    set(gca,'LineWidth',1)

    sub_cnt = sub_cnt + 3;
    
    clear x* y* z*
end



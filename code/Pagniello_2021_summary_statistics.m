%% Summary Statistics

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

%% Modeled Replicas Fields

load('F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\error_type0_model_SCOT.mat','DELTA*','x','y','z','xx','yy','zz');

xy = sqrt(x.^2 + y.^2); xyz = sqrt(x.^2 + y.^2 + z.^2);
xxyy = sqrt(xx.^2 + yy.^2); xxyyzz = sqrt(xx.^2 + yy.^2 + zz.^2);

stats_med_mod = [median(DELTA_x); median(DELTA_y); median(DELTA_z); median(DELTA_xy); median(DELTA_xyz)];
stats_mad_mod = [mad(DELTA_x,1); mad(DELTA_y,1); mad(DELTA_z,1); mad(DELTA_xy,1); mad(DELTA_xyz,1)];
stats_min_mod = [min(DELTA_x); min(DELTA_y); min(DELTA_z); min(DELTA_xy); min(DELTA_xyz)];
stats_max_mod = [max(DELTA_x); max(DELTA_y); max(DELTA_z); max(DELTA_xy); max(DELTA_xyz)];
stats_CCC_mod = [getfield(f_CCC([x xx],0.05),'est'); getfield(f_CCC([y yy],0.05),'est'); getfield(f_CCC([z zz],0.05),'est'); ...
    getfield(f_CCC([xy xxyy],0.05),'est'); getfield(f_CCC([xyz xxyyzz],0.05),'est')];

ind = z <= 8;

stats_med_mod_s = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_mod_s = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_mod_s = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_mod_s = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_mod_s = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

ind = z > 8;

stats_med_mod_d = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_mod_d = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_mod_d = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_mod_d = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_mod_d = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

d = sqrt(x.^2 + y.^2);
ind = d <= array.dist(3);

stats_med_mod_in = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_mod_in = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_mod_in = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_mod_in = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_mod_in = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

ind = d > array.dist(3);

stats_med_mod_out = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_mod_out = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_mod_out = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_mod_out = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_mod_out = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

clear DELTA* x y z ind

%% Data Derived Replica Fields

load('F:\C_files\C_Spring_2021\Localization_Playbacks\Localization_Output\error_type0_data_SCOT.mat','DELTA*','x','y','z','xx','yy','zz');

xy = sqrt(x.^2 + y.^2); xyz = sqrt(x.^2 + y.^2 + z.^2);
xxyy = sqrt(xx.^2 + yy.^2); xxyyzz = sqrt(xx.^2 + yy.^2 + zz.^2);

stats_med_dat = [median(DELTA_x); median(DELTA_y); median(DELTA_z); median(DELTA_xy); median(DELTA_xyz)];
stats_mad_dat = [mad(DELTA_x,1); mad(DELTA_y,1); mad(DELTA_z,1); mad(DELTA_xy,1); mad(DELTA_xyz,1)];
stats_min_dat = [min(DELTA_x); min(DELTA_y); min(DELTA_z); min(DELTA_xy); min(DELTA_xyz)];
stats_max_dat = [max(DELTA_x); max(DELTA_y); max(DELTA_z); max(DELTA_xy); max(DELTA_xyz)];
stats_CCC_dat = [getfield(f_CCC([x xx],0.05),'est'); getfield(f_CCC([y yy],0.05),'est'); getfield(f_CCC([z zz],0.05),'est'); ...
    getfield(f_CCC([xy xxyy],0.05),'est'); getfield(f_CCC([xyz xxyyzz],0.05),'est')];

ind = z <= 8;

stats_med_dat_s = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_dat_s = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_dat_s = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_dat_s = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_dat_s = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

ind = z > 8;

stats_med_dat_d = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_dat_d = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_dat_d = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_dat_d = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_dat_d = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

d = sqrt(x.^2 + y.^2);
ind = d <= array.dist(3);

stats_med_dat_in = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_dat_in = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_dat_in = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_dat_in = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_dat_in = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

ind = d > array.dist(3);

stats_med_dat_out = [median(DELTA_x(ind)); median(DELTA_y(ind)); median(DELTA_z(ind)); median(DELTA_xy(ind)); median(DELTA_xyz(ind))];
stats_mad_dat_out = [mad(DELTA_x(ind),1); mad(DELTA_y(ind),1); mad(DELTA_z(ind),1); mad(DELTA_xy(ind),1); mad(DELTA_xyz(ind),1)];
stats_min_dat_out = [min(DELTA_x(ind)); min(DELTA_y(ind)); min(DELTA_z(ind)); min(DELTA_xy(ind)); min(DELTA_xyz(ind))];
stats_max_dat_out = [max(DELTA_x(ind)); max(DELTA_y(ind)); max(DELTA_z(ind)); max(DELTA_xy(ind)); max(DELTA_xyz(ind))];
stats_CCC_dat_out = [getfield(f_CCC([x(ind) xx(ind)],0.05),'est'); getfield(f_CCC([y(ind) yy(ind)],0.05),'est'); getfield(f_CCC([z(ind) zz(ind)],0.05),'est'); ...
    getfield(f_CCC([xy(ind) xxyy(ind)],0.05),'est'); getfield(f_CCC([xyz(ind) xxyyzz(ind)],0.05),'est')];

clear DELTA* x y z ind

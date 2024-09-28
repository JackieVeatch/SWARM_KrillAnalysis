%% read in EK80, assign to surface & subsurface
% Created by Jackie 01May2023
% major edits 29Aug2023 to fix MLD match and add phyto patch match

% Specify the directory path
dir_path = '/Volumes/T7_Shield/jmv208/SWARM_data/EK80/ACRO_surveys/';

% Get information about the files in the directory
files = dir(fullfile(dir_path, '*.xls'));

% Extract the file names from the returned structure
file_names = {files.name};
EK80.Sv_mean = [];
EK80.NASC = [];
EK80.Height_mean = [];
EK80.Depth_mean = [];
EK80.Date_M = [];
EK80.Lat_M = [];
EK80.Lon_M = [];
EK80.Thickness_mean = [];
EK80.Time_M = [];
EK80.Length = [];

% I don't know why the first 15 filenames are just slightly wrong versions
% of the 15 files I need... but here's a quick fix that avoids that
% existential crisis
for i = 16:length(file_names)
    str = strcat(dir_path, file_names{i});
    survey = readtable(str);
    EK80.Sv_mean = vertcat(EK80.Sv_mean, survey.Sv_mean);
    EK80.NASC = vertcat(EK80.NASC, survey.NASC);
    EK80.Height_mean = vertcat(EK80.Height_mean, survey.Height_mean);
    EK80.Depth_mean = vertcat(EK80.Depth_mean, survey.Depth_mean);
    EK80.Date_M = vertcat(EK80.Date_M, survey.Date_M);
    EK80.Lat_M = vertcat(EK80.Lat_M, survey.Lat_M);
    EK80.Lon_M = vertcat(EK80.Lon_M, survey.Lon_M);
    EK80.Thickness_mean = vertcat(EK80.Thickness_mean, survey.Corrected_thickness);
    EK80.Time_M = vertcat(EK80.Time_M, survey.Time_M);
    EK80.Length = vertcat(EK80.Length, survey.Corrected_length);
    
end


%% convert EK80 weird time format into datetime
date = datetime(num2str(EK80.Date_M), 'InputFormat', 'yyyyMMdd');
dt = datenum(date);
dt_hour = dt + EK80.Time_M; % time variable has fraction of day 
EK80.datenum = dt_hour;

%% load in ACRO data and match krill swarms to above and below MLD
    % most of them will probably be below becuase this happens during the
    % day
    
load /Volumes/T7_Shield/jmv208/SWARM_data/acro_data_reprocessed/ACRO_reprocessed.mat

% match krill swarms in space with ACRO surveys --> more accurate than time
% becuase the ACRO was being towed 20 meters behind the EK80

% % it is suggested by MATHWORKS to use a triangulation before 'dsearchn'
% when using a large amount of points... I don't think we have what is
% considered a 'large' amount of points here so I'm not going to use it for
% now... also I can't get it to work so there's that
% T = delaunayn([ACRO.lon(2:end-1), ACRO.lat(2:end-1)]);
% [k_test, dist_test] = dsearchn([ACRO.lon(2:end-1), ACRO.lat(2:end-1)], T, [EK80.Lon_M, EK80.Lat_M]);
% ind_nomatch_test = find(dist> 0.02);

dt = datetime(EK80.datenum, 'ConvertFrom', 'datenum');
dateComponents = dateshift(dt, 'start', 'day');
days = unique(dateComponents);

mld_all = [];

for i = 1:length(days)
    ind_EK80 = find(EK80.datenum > datenum(days(i)) & EK80.datenum < datenum(days(i))+1);
    ind_ACRO = find(ACRO.mtime > datenum(days(i)) & ACRO.mtime < datenum(days(i))+1);
    
    [k,dist] = dsearchn([ACRO.lon(ind_ACRO)', ACRO.lat(ind_ACRO)'], [EK80.Lon_M(ind_EK80), EK80.Lat_M(ind_EK80)]);
    ind_nomatch = find(dist> 1);
    
    mld_day = ACRO.mld(ind_ACRO);
    mld = mld_day(k);
    mld(ind_nomatch) = NaN;
    
    mld_all = [mld_all, mld];
    
end

EK80.mld = mld_all;

for i = 1:length(EK80.mld)
    if EK80.Depth_mean(i) >= EK80.mld(i)
        EK80.surface(i) = 0;
    else
        EK80.surface(i) = 1;
    end
end


%% use the same methodology as above to assign patch/no phytoplankton patch to each krill swarm

load '/Volumes/T7_Shield/jmv208/ACROBAT/patchID_Adelie_matched_extraDays.mat'

patch_all = [];

for i = 1:length(days)
    ind_EK80 = find(EK80.datenum > datenum(days(i)) & EK80.datenum < datenum(days(i))+1);
    ind_patch = find(patchID_Adelie.timestamp > datenum(days(i)) & patchID_Adelie.timestamp < datenum(days(i))+1);
    
    [k,dist] = dsearchn([patchID_Adelie.lon(ind_patch)', patchID_Adelie.lat(ind_patch)'], [EK80.Lon_M(ind_EK80), EK80.Lat_M(ind_EK80)]);
    ind_nomatch = find(dist> 1);
    
    patch_day = patchID_Adelie.patch(ind_patch);
    patch = patch_day(k);
    patch(ind_nomatch) = NaN;
    
    patch_all = [patch_all, patch];
    
end

EK80.patch = patch_all;

%% do the same as above but matching in time rather than space

patch_all = [];
part_all = [];

for i = 1:length(days)
    ind_EK80 = find(EK80.datenum > datenum(days(i)) & EK80.datenum < datenum(days(i))+1);
    ind_patch = find(patchID_Adelie.timestamp > datenum(days(i)) & patchID_Adelie.timestamp < datenum(days(i))+1);
    
    EK80_time = EK80.datenum(ind_EK80);
    patch_binary = patchID_Adelie.patch(ind_patch);
    int_part = patchID_Adelie.mlparticle(ind_patch);
    patch = [];
    part = [];
    
    for j = 1:length(ind_EK80)
        
        timeDifferences = abs(patchID_Adelie.timestamp(ind_patch) - EK80_time(j));
        ind_closest = find(timeDifferences== min(timeDifferences));
        if timeDifferences(ind_closest) > 1/48 % if the time difference is greater that 30 minutes
            patch(j) = NaN;
            part(j) = NaN;
        else
            patch(j) = patch_binary(ind_closest);
            part(j) = int_part(ind_closest);
        end
        
    end
    
    patch_all = [patch_all, patch];
    part_all = [part_all, part];
    
end

EK80.patch_timeMatch = patch_all;
EK80.part_timeMatch = part_all;

EK80.createdWith = 'define_krill_patches_EK80.m';
EK80.lastEdited = '18Oct2023';

%%%%% FOR SOME REASON THE TIME MATCH GIVES A LOT OF NAN VALUES??
    %%% I think this is becuase the patchID is only Adelie and the EK80
    %%% still has the gentoo transect in it?
    %%% for this reason I think I should use the time matched values






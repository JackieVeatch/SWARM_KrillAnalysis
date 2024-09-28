%% match krill swarm to FTLE value
    % adapted from 'resample_ftle.m' to perform same task on krill swarms
    % created by Jackie 10July2023

load 'EK80.mat';
load '/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset/ftle_LR_binary.mat'

% convert EK80 weird time format into datetime
date = datetime(num2str(EK80.Date_M), 'InputFormat', 'yyyyMMdd');
dt = datenum(date);
dt_hour = dt + EK80.Time_M; % time variable has fraction of day 
EK80.datenum = dt_hour;

%% loop through ACRO observations, pair to ftle

start_codar = min(ftle_LR.time);
end_codar = max(ftle_LR.time);

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
EK80.ftle_match = NaN(size(EK80.datenum));
% create grid

ftle_2d = reshape(ftle_LR.ftle, [4400, 1273]);
[X,Y] = meshgrid(ftle_LR.x, ftle_LR.y);
X = reshape(X, [4400,1]);
Y = reshape(Y, [4400,1]);
ftle_coords = [X,Y];

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((EK80.datenum>=datenum(start))&(EK80.datenum<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = EK80.Lat_M(ind);
        lon = EK80.Lon_M(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = ftle_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            EK80.ftle_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end

%% resample RPD

% modified from 'resample_rpd_edges.m'

% gridded RPD from "autocor_timescale_rpd.m" --> values not coords
load '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/RPD_gridded_season.mat'
load '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/RPD_coordinates.mat'

EK80.rpd_match = NaN(size(EK80.datenum));
rpd_2d = reshape(part_dens_gridded, [667,1511]);

load('/Volumes/T7_Shield/jmv208/SWARM_data/CODAR_zeroed.mat');

for i = 1:length(CODAR_zeroed.time)
    day_sum(i) = sum(CODAR_zeroed.u(:,:,i), 'all');
end
ind = find(day_sum ~=0);
start_codar = CODAR_zeroed.dnum(ind(1));
end_codar = CODAR_zeroed.dnum(ind(end))-3;
% this time index is the same one I used to create trajectories and
% calculate rpd with.

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((EK80.datenum>=datenum(start))&(EK80.datenum<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = EK80.Lat_M(ind);
        lon = EK80.Lon_M(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(RPD_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = rpd_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            EK80.rpd_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end



EK80.edited = 'resample_ftle_EK80.m to include ftle values';



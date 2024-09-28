%% match krill swarm to RPD and FTLE value, randomly generated data
    % adapted from 'resample_ftle_rpd_EK80.m' to perform same task on randomly generated krill swarms
    % created by Jackie 13July2023
    % outfitted for baffin

load '/home/jmv208/krill_code/random_resample.mat';


load '/home/jmv208/SWARM_data/RPD_gridded_season.mat'
load '/home/jmv208/SWARM_data/RPD_coordinates.mat'


load('/home/jmv208/SWARM_data/CODAR_zeroed.mat');
load '/home/jmv208/LCS-Tool-master/demo/ocean_dataset/ftle_LR_binary.mat'
%% loop through EK80 observations, pair to ftle

start_codar = min(ftle_LR.time);
end_codar = max(ftle_LR.time);

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
random_resample.ftle_match = NaN(size(random_resample.time));
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
    ind = find((random_resample.time>=datenum(start))&(random_resample.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample.lat(ind);
        lon = random_resample.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = ftle_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample.ftle_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end

%% resample RPD

% modified from 'resample_rpd_edges.m'

random_resample.rpd_match = NaN(size(random_resample.time));
rpd_2d = reshape(part_dens_gridded, [667,1511]);

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
    ind = find((random_resample.time>=datenum(start))&(random_resample.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample.lat(ind);
        lon = random_resample.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(RPD_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = rpd_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample.rpd_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end



random_resample.edited = 'resample_ftle_rpd_backgroun_EK80.m to include ftle and rpd values';

random_resample_matched = random_resample;

save('random_resample_matched' , 'random_resample_matched');

exit
clear; clc; close all;

load('RCVRT_data.mat')

% --- Params
c = physconst('LightSpeed'); % speed of light

% pre-allocate cell size
svData = cell(length( RCVRT_data.gpsTime(:,1)),1);

for i = 1:length(RCVRT_data.gpsTime(:,1))
    
    % find sats in view
    [~, satsInView] = find(~isnan(RCVRT_data.obs.L1.psr(i,2:end)));
    
    satsInView = satsInView + 1;
    
    % extract psr dopp time
    svData{i}.gpsTime = RCVRT_data.gpsTime(i,1);
    svData{i}.psr(:,1) = RCVRT_data.obs.L1.psr(i,satsInView);
    svData{i}.dopp(:,1) = RCVRT_data.obs.L1.dopp(i,satsInView);
    svData{i}.satsInView = satsInView; % store sats
    % calc pos and vel
    transitTime = svData{i}.psr/c;
    transmitTime = RCVRT_data.gpsTime(i) - transitTime;
    for j = 1:length(transitTime)
        [svData{i}.pos(j,:), svData{i}.vel(j,:), svData{i}.clkCorr(j,1)] = calc_sv_pos2(RCVRT_data.ephem(satsInView(j),:), transmitTime(j), transitTime(j));
    end 
    
end 

save('satData.mat', 'svData')
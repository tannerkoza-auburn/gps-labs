clear
clc
close all

dir = fileparts(which(mfilename)); % add all subfolders at run
addpath(genpath(dir));

% truth given from surveyed position
ecef_t = [423203.359; -5361678.541; 3417280.681];

%% Initialization & Data Import

load('satData.mat')


k = length(svData); % number of time steps
ts = 1:k; % time step vector

estPos = zeros(3, k); % preallocate position estimates
estVel = zeros(3, k); % preallocate velocity estimates

estClockBias = 0; % initial clock bias estimate (m)

measSig = 0.5; % receiver accuracy (m)

wgs84 = wgs84Ellipsoid('meter'); % reference spheroid

%% State Estimate

% --- rotation ECEF to ENU for DOP rotation
lat0 = 32.606460; % toomer's lat
lon0 = -85.481854; % toomer's lon
h0 = 214; % m
R = [-sind(lon0) -sind(lat0)*cosd(lon0) cosd(lat0)*cosd(lon0);
    cosd(lon0)      -sind(lat0)*sind(lon0)  cosd(lat0)*sind(lon0);
        0           cosd(lat0)          sind(lon0)];

rcvr = gnssReceiver;

[svAz, svEl] = deal( zeros(32,length(ts)) );
satsInViewIdx = zeros(32,length(ts));
for i = ts
    
    % --- Unpack Data
    svPos = svData{i}.pos;
    svVel = svData{i}.vel;
    svPsr = svData{i}.psr;
    svDopp = svData{i}.dopp;
    satsInView = svData{i}.satsInView; % PRNs of sats in view
    satsInViewIdx(svData{i}.satsInView,i) = 1;
    clkCorr = svData{i}.clkCorr;
    time(i) = svData{i}.gpsTime;
   
    stateEst = rcvr.ddp3D(svPsr, svPos, svVel, estPos(:,i), estClockBias, measSig);
    
    % correct state estimate
    estPos(:,i) = stateEst.Pos;
    estVel(:,i) = stateEst.Vel;
    
    if i < ts
    % initialize for next time step
    estPos(:,i+1) = stateEst.Pos;
    estVel(:,i+1) = stateEst.Vel;
    end 
    
    % calculate az/el positions at this time step:
    for j = 1:length(svPsr)        
        lla = ecef2lla(estPos(:,1)');
        lat_aer = lla(1);
        lon_aer = lla(2);
        h_aer = lla(3);
        % 3D matrix of az and el of satellites at time, t
        [svAz(satsInView,i), svEl(satsInView,i), ~] = ecef2aer(svPos(:,1), svPos(:,2), svPos(:,3), lat_aer, lon_aer, h_aer, wgs84);        
    end 
    
    DOP_ENU = R'*stateEst.DOP(1:3,1:3)*R;
    DOP_E(i) = sqrt(DOP_ENU(1,1));
    DOP_N(i) = sqrt(DOP_ENU(2,2));
    DOP_U(i) = sqrt(DOP_ENU(3,3));

end

time = time - time(1); % zero out time

% % --- true position for comparison
% [E_t, N_t, U_t] = ecef2enu(ecef_t(1), ecef_t(2), ecef_t(3), lat0, lon0, h0, wgs84);
% 
% % --- Part A
% % find where sats come in and out and create diff figures
% idx_1 = find(satsInViewIdx(1,:)); idx_1 = idx_1(1); % sat 1 comes into view
% idx_14 = find(~satsInViewIdx(14,:)); idx_14 = idx_14(1); % sat 14 leaves view
% 
% idx_fig1 = min([idx_1 idx_14]); % first data set from 1:idx_fig1
% idx_fig2 = max([idx_1 idx_14]); % second data set from idx_fig1:idx_fig2
% 
% figure()
% skyPlot(svAz(find(svAz(:,1)),1:idx_fig1-1), svEl(find(svAz(:,1)),1:idx_fig1-1), find(satsInViewIdx(:,1)));
% title(sprintf('Satellites in View: 0 to %.2f sec', time(idx_fig1-1)));
% 
% figure()
% skyPlot(svAz(find(svAz(:,idx_fig1)),idx_fig1:idx_fig2-1), svEl(find(svAz(:,idx_fig1)),idx_fig1:idx_fig2-1), find(satsInViewIdx(:,idx_fig2-1)));
% title(sprintf('Satellites in View: %.f to %.2f sec', time(idx_fig1), time(idx_fig2-1)));
% 
% figure()
% skyPlot(svAz(find(svAz(:,idx_fig2)),idx_fig2:end), svEl(find(svAz(:,idx_fig2)),idx_fig2:end), find(satsInViewIdx(:,end)));
% title(sprintf('Satellites in View: %.f to %.2f sec', time(idx_fig2), time(end)));




% % --- Part B
% figure
% lla_t = ecef2lla(ecef_t');
% lla = ecef2lla(estPos');
% geoplot(lla_t(1), lla_t(2), '*b')
% hold on
% geoplot(lla(1,1), lla(1,2),'*r')
% legend('Provided Initial Positon', 'Estimatated Initial Position')
% geobasemap satellite
% 
% 
% % --- Part D
% [E, N, U] = ecef2enu(estPos(1,:), estPos(2,:), estPos(3,:), lat0, lon0, h0, wgs84);
% err = [E - E_t; N - N_t; U - U_t];
% figure()
% subplot(3,1,1)
% plot(time,err(1,:));
% title('East Error')
% ylabel('East Error (m)')
% xlabel('Time (s)')
% subplot(3,1,2)
% plot(time,err(2,:));
% title('North Error')
% ylabel('North Error(m)')
% xlabel('Time (s)')
% subplot(3,1,3)
% plot(time,err(3,:));
% title('Up Error')
% ylabel('Up Error (m)')
% xlabel('Time (s)')
% 
% figure()
% plot(E, N)
% hold on
% xlabel('East (m)')
% ylabel('North (m)')
% plot(0, 0, '*r', 'MarkerSize', 20);
% legend('Estimated Position', 'Toomers Corner')
% 
% 
% figure()
% subplot(3,1,1)
% plot(time,DOP_E, 'linewidth', 2);
% title('East DOP')
% ylabel('East DOP')
% xlabel('Time (s)')
% subplot(3,1,2)
% plot(time,DOP_N, 'linewidth', 2);
% title('North DOP')
% ylabel('North DOP')
% xlabel('Time (s)')
% subplot(3,1,3)
% plot(time,DOP_U, 'linewidth', 2);
% title('Up DOP')
% ylabel('Up DOP')
% xlabel('Time (s)')
% 
% 
% figure()




[velEast,velNorth,velUp] = ecef2enuv(estVel(1,:),estVel(2,:),estVel(3,:),lat0,lon0);


figure
subplot(3,1,1)
plot(ts,velEast)
subplot(3,1,2)
plot(ts,velNorth)
subplot(3,1,3)
plot(ts,velUp)




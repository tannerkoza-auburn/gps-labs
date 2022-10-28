clear
clc
close all

dir = fileparts(which(mfilename)); % add all subfolders at run
addpath(genpath(dir));

%% Initialization & Data Import

load('satDataP2.mat')

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
R = [-sind(lon0) -sind(lat0)*cosd(lon0) cosd(lat0)*cosd(lon0);
    cosd(lon0)      -sind(lat0)*sind(lon0)  cosd(lat0)*sind(lon0);
        0           cosd(lat0)          sind(lon0)];



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
    DOPP(:,i) = svDopp;
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
        lat0 = lla(1);
        lon0 = lla(2);
        h0 = lla(3);
        % 3D matrix of az and el of satellites at time, t
        [svAz(satsInView,i), svEl(satsInView,i), ~] = ecef2aer(svPos(:,1), svPos(:,2), svPos(:,3), lat0, lon0, h0, wgs84);        
    end 
    
    DOP_ENU = R'*stateEst.DOP(1:3,1:3)*R;
    DOP_E(i) = sqrt(DOP_ENU(1,1));
    DOP_N(i) = sqrt(DOP_ENU(2,2));
    DOP_U(i) = sqrt(DOP_ENU(3,3));

    
    clkBias(i) = stateEst.ClockBias;
    clkDrift(i) = stateEst.ClockDrift;
    
end

time = time - time(1); % zero out time



% --- Part D
lat0 = 32.606460; % toomer's lat
lon0 = -85.481854; % toomer's lon
h0 = 214; % height of auburn in km
[E, N, U] = ecef2enu(estPos(1,:), estPos(2,:), estPos(3,:), lat0, lon0, h0, wgs84);
figure()
subplot(3,1,1)
plot(time,E);
title('East')
ylabel('East (m)')
xlabel('Time (s)')
subplot(3,1,2)
plot(time,N);
title('North')
ylabel('North (m)')
xlabel('Time (s)')
subplot(3,1,3)
plot(time,U);
title('Up')
ylabel('Up (m)')
xlabel('Time (s)')

plot(E, N)
hold on
xlabel('East (m)')
ylabel('North (m)')
plot(0, 0, '*r', 'MarkerSize', 20);
legend('Estimated Position', 'Toomers Corner')


figure()
subplot(3,1,1)
plot(time,DOP_E, 'linewidth', 2);
title('East DOP')
ylabel('East DOP')
xlabel('Time (s)')
subplot(3,1,2)
plot(time,DOP_N, 'linewidth', 2);
title('North DOP')
ylabel('North DOP')
xlabel('Time (s)')
subplot(3,1,3)
plot(time,DOP_U, 'linewidth', 2);
title('Up DOP')
ylabel('Up DOP')
xlabel('Time (s)')

sig_psr = 0.5; % (m)
muE = mean(sig_psr*DOP_E);
muN = mean(sig_psr*DOP_N);
muU = mean(sig_psr*DOP_U);
sigE = std(sig_psr*DOP_E);
sigN = std(sig_psr*DOP_N);
sigU = std(sig_psr*DOP_U);
figure()
subplot(3,1,1)
plot(time,sig_psr*DOP_E, 'linewidth', 2);
hold on
title('East DOP Estimated Error')
ylabel('East Err (m)')
xlabel('Time (s)')
subplot(3,1,2)
plot(time,sig_psr*DOP_N, 'linewidth', 2);
title('North DOP Estimated Error')
ylabel('North Err (m)')
xlabel('Time (s)')
subplot(3,1,3)
plot(time,sig_psr*DOP_U, 'linewidth', 2);
title('Up DOP Estimated Error')
ylabel('Up Err (m)')
xlabel('Time (s)')

figure()
subplot(2,1,1)
plot(time, clkBias, 'linewidth', 2)
subplot(2,1,2)
plot(time, clkDrift, 'linewidth', 2)


lla = ecef2lla(estPos');
[velEast,velNorth,velUp] = ecef2enuv(estVel(1,:),estVel(2,:),estVel(3,:),lla(2,1),lla(2,2));

speed = sqrt(velEast.^2 + velNorth.^2 + velUp.^2);

course = 0;

for i = 2:length(velEast)
    if speed(i) > 0.5
        course(i) = atan2d(velEast(i), velNorth(i));
   
    else
        course(i) = course(i-1);
     
    end
end

figure
geoplot(lla(2:end,1), lla(2:end,2),'r.')
geobasemap satellite

figure
subplot(3,1,1)
plot(ts,velEast)
subplot(3,1,2)
plot(ts,velNorth)
subplot(3,1,3)
plot(ts,velUp)




function stateEst = gnssPosition(rho, rhodot, svClockCorr, svPos, svVel, estPos, estClockBias, rcvrSigma)
% DESCRIPTION: This function produces a global position estimate using 
% least squares error estimation & the Newton-Raphson provided satellite
% positions and pseudoranges. The number of least squares iterations
% used is also provided.
% PARAMS:
%       rho: column vector of pseudoranges
%       svPos: nxm matrix of satellite(s) positions
% OUTPUT:
%       A struct full of the following fields:
%       estPos: estimated position solution
%       estClockBias: estimated clock bias
%       P: estimated covaraince of our position error (dPos) estimate
%       itr: number of iterations taken to converge
% AUTHOR: Tanner Koza

%% Initialization
    
    % Handle Input Dimensions
    numMeas = length(rho);
    sz_rho = size(rho);
    sz_rhodot = size(rhodot);
    sz_svPos = size(svPos);
    sz_svVel = size(svVel);

    if sz_rho(1) ~= numMeas

        rho = rho';
    end

    if sz_rhodot(1) ~= numMeas

        rhodot = rhodot';
    end

    if sz_svPos(2) ~= numMeas

        svPos = svPos';
    end

    if sz_svVel(2) ~= numMeas

        svVel = svVel';

    end


    % Initialize Position & Clock Bias Estimate
    estX = [estPos; estClockBias; zeros(4,1)];

    % Initialize Iteration Count
    itr = 0;

    % Signal Velocity
    sig_freq = 1575.42e6; % carrier signal frequency

    % Initialize Convergence Criterion
    conv_th = 10e-8;

    % Calculate Measurement Variance
    var = rcvrSigma^2;

    % Convert rhodot to m/s
    C = physconst('LightSpeed');

    rhodot = rhodot * -(C/sig_freq);

    % Convert SV Clock Corrections to Meters
    svClockCorr = svClockCorr * C;

%% Estimation

    % Least Squares & Newton-Raphson
    while true

        unitVecs = gnssUnitVector(svPos, estPos);
    
        y = gnssMeasVector(rho, rhodot, svClockCorr, unitVecs, svPos, svVel, estPos, estClockBias);

        G = gnssGeomMatrix(unitVecs);

        dX = (G' * G)^-1 * G' * y;

        estX = estX + dX;

        estPos = estX(1:3);
        estVel = dX(5:7);
        estClockBias = estX(4);
        estClockDrift = dX(8);

        itr = itr + 1;

        if  norm(dX(1:4)) <= conv_th
            break
        end

    end  

    % Predicted Error Estimate Covariance Matrix
    DOP = ( G' * G )^-1;
    P = var * DOP;

    % Solution Structure Population
    stateEst.Pos = estPos;
    stateEst.Vel = estVel;
    stateEst.ClockBias = estClockBias;
    stateEst.ClockDrift = estClockDrift;
    stateEst.DOP = DOP;
    stateEst.P = P;
    stateEst.itr = itr;

end
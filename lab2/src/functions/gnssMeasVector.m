function y = gnssMeasVector(rho, rhodot, svClockCorr, unitVecs, svPos, svVel, estPos, estClockBias)
% DESCRIPTION: This function produces the measurement vector used in the
% least squares estimation of global position using satellite pseudoranges.
% PARAMS:
%       rho: column vector of pseudoranges
%       svPos: nxm matrix of satellite(s) positions
%       estPos: column vector of current estimated user position
%       estClockBias: scalar clock bias (m)
% OUTPUT:
%       y: GNSS measurement vector
% NOTES: 
%   - estClockBias is range bias that stems from the user's clock
%   bias. ( Ti - Tu ) is the clock offset between the user's and
%   satellite's clocks.
% AUTHOR: Tanner Koza, M.E. (Master of Engineering) Candidate

%% Initialization

    % Preallocate Estimated Pseudorange Vector 
    numMeas = length(rho);
    rhohat = zeros(numMeas,1);
    rhodothat = zeros(numMeas,1);

%% Measurement Vector Population

    % Calculate Estimated Pseudorange
    for i = 1:numMeas

        rhohat(i) = sqrt( ( svPos(1,i) - estPos(1) )^2 ...
            + ( svPos(2,i) - estPos(2) )^2 ...
            + ( svPos(3,i) - estPos(3) )^2) + estClockBias - svClockCorr(i);

        rhodothat(i) = unitVecs(i,1) * svVel(1,i) + ...
                       unitVecs(i,2) * svVel(2,i) + unitVecs(i,3) * svVel(3,i);

    end

    y = [rho; rhodot] - [rhohat; rhodothat];

end
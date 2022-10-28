function  U = gnssUnitVector(svPos, estPos)
% DESCRIPTION: This function produces the unit vectors 
% used in constructing the measurement vector and geometry matrix 
% in the least squares estimation of global position using
% satellite pseudoranges.
% PARAMS:
%       svPos: nxm matrix of satellite(s) positions
%       estPos: column vector of current estimated user position
% OUTPUT:
%       U: current position unit vectors to satellites
% AUTHOR: Tanner Koza


%% Initialization

    % Preallocate Satellite Unit Vectors 
    numMeas = length(svPos);

    uhat_x = zeros(numMeas,1);
    uhat_y = zeros(numMeas,1);
    uhat_z = zeros(numMeas,1);

%% Unit Vector Calculation

     % Calculate Satellite Unit Vectors
    for i = 1:numMeas

        r = sqrt( ( svPos(1,i) - estPos(1) )^2 ...
            + ( svPos(2,i) - estPos(2) )^2 ...
            + ( svPos(3,i) - estPos(3) )^2);

        uhat_x(i) = ( svPos(1,i) - estPos(1) )/ r;

        uhat_y(i) = ( svPos(2,i) - estPos(2) )/ r;

        uhat_z(i) = ( svPos(3,i) - estPos(3) )/ r;

    end

    U = [uhat_x uhat_y uhat_z];

end
function  G = gnssGeomMatrix(unitVecs)
% DESCRIPTION: This function produces the geometry/observation matrix 
% used in the least squares estimation of global position using
% satellite pseudoranges.
% PARAMS:
%       unitVecs: mx3 matrix of unit vectors to satellite(s) positions
% OUTPUT:
%       G: GNSS geometry matrix
% AUTHOR: Tanner Koza

%% Initialization

    % Preallocate Satellite Unit Vectors 
    numMeas = length(unitVecs);

    uhat_x = unitVecs(:,1);
    uhat_y = unitVecs(:,2);
    uhat_z = unitVecs(:,3);

%% Unit Vector Calculation

    G = [-uhat_x -uhat_y -uhat_z ones(numMeas,1) zeros(numMeas,4);
        zeros(numMeas,4) -uhat_x -uhat_y -uhat_z ones(numMeas,1)];

end
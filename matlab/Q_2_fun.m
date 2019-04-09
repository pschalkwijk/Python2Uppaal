function Q = Q_2_fun(n,m)
% This function generates an isotropic covering for 2-D case. Each conic
% region is described by { x | x^T Qs x >= 0 }. The number of conic regions
% covering the half space is m. This function needs two input variables n
% and m. The input variable n represents dimension of the state space. The
% input variable m represents the number of angular subdivisions of [0,pi],
% i.e. only w.r.t. the top half state-space. 
% NOTE: In The 0 angle corresponds to the x positive (in the (x,y)-plane)
% and increases in the counterclockwise direction.
% The output variable Q is a one-column cell array, where each element
% contains a 2x2 symmetric matrix Qs.

q = m;               % number of regions covering the half space

%%%%%%%%%%%%%%%%%%%%% Calculating all the regions bounds in 2D
thet_region = (1:m);
thet_min = (thet_region-ones(size(thet_region)))*pi/m;
thet_max = thet_region*pi/m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the output variable
Q = cell(q,1);

% create matrix Qs for all s
for j = 1:q

    clear a_1 a_2    % Note: a_1^T x>=0 and a_2^T x>=0

    a_1 = [-sin(thet_min(1,j)) cos(thet_min(1,j))]';
    a_2 = [sin(thet_max(1,j)) -cos(thet_max(1,j))]';
    Q{j,1} = a_1*a_2'+a_2*a_1';

end

end
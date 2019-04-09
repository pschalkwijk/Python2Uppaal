%% Define nD Region Polyhedra
clear 
clc
close all
%% Add Toolbox Directories

% Add YALMIP Directories
YALMIP_dir = strcat(cd,'\YALMIP');
addpath(genpath(YALMIP_dir));

%% Parameters
n = 3;
m = 4;
q = 2*m^(n-1);

%% Initialize regions for state-space abstraction
Q = Q_2_fun(2,m); % Generate Q matrices for two-dimensional planes, given the number of subdivions m

% Obtain nD regions by combining 2D regions (projections of nD regions on (n-1) 2D planes)
Two_dim_div = cell((n-1),1);

for i = 1:(n-2)
    Two_dim_div{i,1} = [1:m];                   % ((n-2) x m - 2D regions
end

Two_dim_div{n-1,1} = [1:2*m];                   % 2*m - 2D regions

AllRegions = combvec(Two_dim_div{1:(n-1),1});   % All possible combinations of 2D regions (each combination corresponding to an nD region

%% Create Polyhedra

Polyhedra = cell(q,1);
medpoints = cell(q,1);

Beams = cell((n-1),1);
r = 100;
other_coor = [r -r];
V = 3*2^(n-2);          % Number of vertices per beam

    for c = 1:(n-2)
        other_coor = combvec(other_coor,[r -r]); % Generate all possible combinations of (n-2) coordinate endpoints
    end
%%
    
for k = 1:(q/2) % Loop on regions
    
    region_vec = AllRegions(:,k);   % Pick 2D projection indices of current region
    
    for b = 1:(n-1)                 % Loop on projections
        proj = region_vec(b);
        
        % Find min. and max. angle for the current projection (in it's own
        % 2D plane)
        theta_min = (region_vec(b,1) - 1) * (pi/m);
        theta_max =  region_vec(b,1) * (pi/m);
        
        % Calculate projection coordinates
        if b == 1
            x1min = cos(theta_min)
            x2min = sin(theta_min)
            x1max = cos(theta_max)
            x2max = sin(theta_max)
        else
            x1min = cos(theta_min-pi/2)
            x2min = sin(theta_min-pi/2)
            x1max = cos(theta_max-pi/2)
            x2max = sin(theta_max-pi/2)
        end
        
        % Write beam coordinates
        coordinate = zeroes(1,n)
            for c = 1:(b-1)
    end
    
 
        
    end

    
end

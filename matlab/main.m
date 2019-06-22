%% General Parameters
close all;
clear all;
clc;

%tic
%% Add Toolbox Directories

% Add YALMIP Directories
%YALMIP_dir = strcat(cd,'\YALMIP');
%addpath(genpath(YALMIP_dir));

% Add MPT Directory
%addpath('H:\Documents\MPT');                           
 

%% Start parallel pool
gcp();

%% Define Control System for Abstraction
% Control loop 1 from Dieky

A = [0 1; 2 -3];
B = [0;1];
K = [1, -4];
% Closed-loop: x_dot = Ax + Bu = Ax - BKx_k
%% Abstraction Parameters
n = length(A(:,1));   % State space dimension

% Considered design parameters
N_conv = 5;           % Order of Taylor series approximation of Phi
sigma_bar = 2;        % Upper limit for sampling time global lower bound (used in step 1)
l = 100;              % Number of subdivisions in the interval [0,sigma_bar]

% m must be an even number
m = 4;               % Number of considered subdivisions for the angle(s) theta (ranging from 0 to pi), must be an even number! 
q = m^(n-1);         % Number of regions (for half the state space)

alpha = 0.05;        % Triggering co�fficient

%% Initialize regions for state-space abstraction
global Q
Q = Q_2_fun(2,m); % Generate Q matrices for two-dimensional planes, given the number of subdivions m

% Obtain nD regions by combining 2D regions (projections of nD regions on (n-1) 2D planes)
Two_dim_div = cell((n-1),1);

for i = 1:(n-2)
    Two_dim_div{i,1} = [1:m];                   % ((n-2) x m - 2D regions
end

Two_dim_div{n-1,1} = [1:2*m];                   % 2*m - 2D regions

global AllRegions
AllRegions = combvec(Two_dim_div{1:(n-1),1});  % All possible combinations of 2D regions (each combination corresponding to an nD region
      
%% Time steps and bounds for each step

% Step 1
tau_min = 0;        % Time to start line search on tau_opt (lower bound on inter-sample time for the entrie state space, without nu) from                   
del_tau = 0.001;    % Time step size in the line search on tau_opt

% Step 2
del_tau_2 = 0.001;  % Time step size in the line search on nu tau_opt_nu (lower bound on inter-sample time for the entrie state space, with nu)
del_sig = 0.001;    % sigma_prime increment

% Step 3
del_tau_3 = 0.001;  % Time step size in the line search on the lower inter-sample time bounds
global epsilon_tol sedumi_eps
epsilon_tol = 0.1;
% Tolerance for constraint on epsilon's for solving LMI's. 
% Should be larger than zero. Setting this to small could cause numerical 
% issues in thez LMI's epsilons, always check if these are positive 
% (for both upper and lower time bound calculations) after the abstraction 
% is completed!!!
sedumi_eps = 1e-5; 
% Precision for sedumi solver; default = 1e-5. 
% Set this precision smaller than epsilon_tol to prevent numerical issues 
% in the LMI's

% Step 4
global sigma_max             % Upper bound to sampling time upper (and lower) bounds
l_star = 2/0.01;              % Number of subdivisions in the interval [0,sigma_max]
del_tau_4 = 0.001;          % Time step size in the line search on the upper inter-sample time bounds

% Simulation          % Initial condition for the simulation
time_end = 20;              % Simulation end time
ts = 0.001;                 % Simulation time step
x_0 = [1; 100];

%% Run abstraction and sampling time bound calculation steps

% Q-MATRICES METHOD
% Make abstraction and calculate sampling time bounds using the Q-matrices
tic
step1_etc
global nu
ctime = toc;
fprintf('finished step 1 after %02.0f:%02.0f:%02.0f\n', [floor(ctime/3600), floor(mod(ctime, 3600)/60), round(mod(ctime, 60))]);
step2_Nu_etc
ctime = toc;
fprintf('finished step 2 after %02.0f:%02.0f:%02.0f\n', [floor(ctime/3600), floor(mod(ctime, 3600)/60), round(mod(ctime, 60))]);
%%
fprintf('tau_opt = %02.3f\n', tau_opt_nu)
Tau_s_opt = step3_nD_Q_etc_lowerbounds(m, sigma_bar, l, N_conv, ...
                                alpha, A, B, K, n, tau_opt_nu, del_tau_3, ...
                                AllRegions, Q, nu, sedumi_eps, epsilon_tol);
ctime = toc;
fprintf('finished lower bounds after %02.0f:%02.0f:%02.0f\n', [floor(ctime/3600), floor(mod(ctime, 3600)/60), round(mod(ctime, 60))]);
%%
sigma_max = 2; 
Tau_s_max = step4_nD_Q_etc_upperbounds(m, l_star, N_conv, alpha, A, B, K, n, Tau_s_opt, del_tau_4, ...
                                         AllRegions, Q, sigma_max, nu, sedumi_eps, epsilon_tol);
% Give time that the inter-sample time bound calculations took to complete in seconds
ctime = toc;
% Give time that the inter-sample time bound calculations took to complete in seconds
fprintf('finished upper bounds after %02.0f:%02.0f:%02.0f\n', [floor(ctime/3600), floor(mod(ctime, 3600)/60), round(mod(ctime, 60))]);



%% Perform reachability analysis
% (Only works if Q-matrices method is used!!!)
% Also: not working for 4D or higher-D yet

[Reachable_regions_regDetQ, Reachable_regions] = step5_nD_Reachability_MPT(m, n, AllRegions, Tau_s_opt, Tau_s_max, A, B, K);

reachability_calculation_time = toc;
% reachability_calculation_time_hrs = reachability_calculation_time/60
 ctime = toc;
% Give time that the inter-sample time bound calculations took to complete in seconds
fprintf('finished reachability analysis after %02.0f:%02.0f:%02.0f\n', [floor(ctime/3600), floor(mod(ctime, 3600)/60), round(mod(ctime, 60))]);

plot_transition_table

%% Run simulation
fig_Sim_Q_etc

%% Plot results
plot_figures

%% Save data
save("abstraction")



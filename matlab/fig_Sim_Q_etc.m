%% Simulation of the systems based on the derived sampling function

% Simulation
% x_0 = [6.21; -0.2; 2.37];
% time_end = 30;
% ts = 0.001;

size_x = ceil(time_end/ts);              % size of results based on time span
x_resp = zeros(n,size_x);                % states evolution, resp for response
x_norm_resp = zeros(1,size_x);           % norm of states
e_resp = zeros(n,size_x);                % errors of states
e_norm_resp = zeros(1,size_x);           % norm of errors
u_0 = zeros(1,size_x);                   % input

t_reset = 0;                             % resetter of time for the evolution of states
t_period = -ts;                          % inner loop period (inter-sample time)   
x_0_mid = x_0;                           % inner loop initial state (sampled state)

trig_time = [];                          % triggering time sequence
i_trig =[];                              % triggering index sequence
x_trig =[];                              % triggering state sequence

reg_Seq = reg_Det_Q(x_0,m);              % region sequence at triggering times 
%reg_Det_Q calculates the region that x_0 lies in, m is the number of partitions for each angular coordinate.  
                                                                                    

% Numerical integration
% Able to handle singular matrices
    eAt = @(t) expm(A*t);               
    
    int_eAx = @(x) integral(eAt,0,x,'ArrayValued',true);
    
clear k t
k=1;                                     % index for states, reset when the state is sampled

for t=1:size_x
    t_period = t_period+ts;         % the time elapsed since the last triggering
    x_k_1 = (eye(n)+int_eAx(((t-1)*ts-t_reset))*(A-B*K))*x_0_mid; % the current state
    x_resp(:,k) = x_k_1;            % collect the current state
    x_norm_resp(1,k) = norm(x_k_1); % collect norm of the current state
    e_k = x_0_mid-x_k_1;            % compute current error = difference between ''current state'' and ''sampled state''
    e_resp(:,k) = e_k;              % collect current error
    e_norm_resp(1,k) = norm(e_k);   % collect norm of current error
    u_0(k) = -K*x_0_mid;          % input
    if norm(e_k)>alpha^(0.5)*norm(x_k_1)  % triggering mechanism
        x_0_mid = x_k_1;                  % initial state reset --> because of triggering (x_0_mid is the sampled state)
        trig_time = [trig_time t_period]; % collect the inter-sample time
        reg_Seq = [reg_Seq reg_Det_Q(x_0_mid,m)]; % collect the region of the sampled state
        t_period = 0;                    % inner loop time reset --> because of triggering
        t_reset = (t-1)*ts;              % backward time shift --> because of triggering
        i_trig = [i_trig, t];            % collect triggering index    
        x_trig = [x_trig, x_0_mid];      % collect state at triggering (sampled state)
    end
    k = k+1; 
end

%% Finding the tau_min and tau_max corresponding to the region of the states at triggering time

Tau_min_sec = [Tau_s_opt Tau_s_opt];
Tau_max_sec = [Tau_s_max Tau_s_max];

clear Min_reg_seq Max_reg_seq
yy = size(trig_time,2);
Min_reg_seq = zeros(1,yy);
Max_reg_seq = zeros(1,yy);
for j=1:yy
    Min_reg_seq(1,j) = Tau_min_sec(1,reg_Seq(1,j));
    Max_reg_seq(1,j) = Tau_max_sec(1,reg_Seq(1,j)); 
end

% NOTE: reg_Seq contains the region for x_0, which is not counted as a
% triggering. The bounds associated to the region of x_0 should hold for
% the first triggering that occurs, since for each triggering the bounds of
% the PREVIOUS sampled state (at the previous triggering) should hold.





 

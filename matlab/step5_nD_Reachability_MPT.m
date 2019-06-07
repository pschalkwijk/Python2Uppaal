%% Reachability Analysis nD
%  This script performs a reachability anaylis for an nD system following
%  the approach described in 'Computing Polyhedral Approximations to Flow
%  Pipes for Dynamic Systems' by Chutinan and Krogh.
function [Reachable_regions_regDetQ, Reachable_regions] = ...
    step5_nD_Reachability_MPT(m, n, AllRegions, Tau_s_opt, Tau_s_max, A, B, K)
%% Parameters
%n = length(A(:,1));         % State-space dimension
%m = 10;                     % Number of divisions for half of the two-dimensional state-spaces
q = 2*m^(n-1);

%% Obtain nD regions by combining 2D regions (projections of nD regions on (n-1) 2D planes)
%Two_dim_div = cell((n-1),1);

%for i = 1:(n-2)
%    Two_dim_div{i,1} = [1:m];                   % ((n-2) x m - 2D regions
%end

%Two_dim_div{n-1,1} = [1:2*m];                   % 2*m - 2D regions

%AllRegions = combvec(Two_dim_div{1:(n-1),1});   % All possible combinations of 2D regions (each combination corresponding to an nD region)
%% Create Region Polyhedra (3D)

Polyhedra = cell(q,1);
medpoints = cell(q,1);

figure
hold on
for i = 1:(q/2) % Loop on regions 
    
    region2D = AllRegions(:,i); % Contains indices for Q-matrices corresponding to the i-th region
    
    theta_min = zeros((n-1),1);
    theta_max = zeros((n-1),1);
    theta_med = zeros((n-1),1);
    
    for j = 1:(n-1) % Loop on projections onto 2D (x_i,x_{i+1})-planes
        % Find min. and max. angle for each 2D projection
        % NOTE: these all lie within [0,pi]!!!
        theta_min(j,1) = (region2D(j,1) - 1) * (pi/m);
        theta_max(j,1) =  region2D(j,1) * (pi/m);
        
        % Find average of each min. and max. angle
        theta_med(j,1) = (theta_min(j,1)+theta_max(j,1))/2;
        
    end
    
    theta_minmax = cell((n-1),1);
    
    for k = 1:(n-1) % Loop on projections onto 2D (x_i,x_{i+1})-planes
        % Store min. and max. angle for each projection
        theta_minmax{k,1} = [theta_min(k) theta_max(k)];      
    end
    
    % Take all possible combinations of min. and max. angles
    angle_combs = combvec(theta_minmax{1:(n-1),1});
    
    V = zeros(2^(n-1),n); % To store the region's 2^(n-1) vertices 
    h_end = 2^(n-1);
    for h = 1:h_end  % Loop on vertices for polyhedron
       % Calculate the region vertices
       % NOTE: for each angle, except theta1, (pi/2) is substracted
       % so that theta_1 lies within [0,pi] and theta_i lies within [-pi/2,pi/2]
       switch n
           case 2
                    V(h,1) = cos(angle_combs(1,h));
                    V(h,2) = sin(angle_combs(1,h));
           case 3
                    if abs(angle_combs(2,h)-(pi/2)) ~= (pi/2)
                        V(h,1) = cos(angle_combs(1,h));
                        V(h,2) = sin(angle_combs(1,h));
                        V(h,3) = sin(angle_combs(1,h))*tan(angle_combs(2,h)-(pi/2));
                    else
                        V(h,1) = 0;
                        V(h,2) = 0;
                        V(h,3) = sign(angle_combs(2,h)-(pi/2));                  
                    end
                    
           otherwise  
                     V(h,1) = cos(angle_combs(1,h));
                     V(h,2) = sin(angle_combs(1,h));
                     
                     for c = 2:(n-1)
                         if abs(angle_combs(c,h)-(pi/2)) ~= (pi/2)
                            V(h,c+1) = V(h,c)*tan(angle_combs(c,h)-(pi/2));
                         else
                             for cc = 1:c                            
                            V(h,cc) = 0; 
                             end
                             
                            V(h,c+1) = sign(angle_combs(c,h)-(pi/2));
                            
                         end
                     end
       end
              
    end
    
    V = unique(V,'rows');   % Use only unique vertices
    
    x_med = zeros(1,n);     % To store middle point of each region
    
    % Calculate middle point of each region
    % NOTE: for each angle, except theta1, (pi/2) is substracted
    % so that theta_1 lies within [0,pi] and theta_i lies within [-pi/2,pi/2]
    switch n
        case 2
                x_med(1,1) = cos(theta_med(1));
                x_med(1,2) = sin(theta_med(1));
%         case 3
%                 if abs(theta_med(2)-(pi/2)) ~= (pi/2)
%                     x_med(1,1) = cos(theta_med(1));
%                     x_med(1,2) = sin(theta_med(1));
%                     x_med(1,3) = sin(theta_med(1))*tan(theta_med(2));
%                 else
%                     x_med(1,1) = 0;
%                     x_med(1,2) = 0;
%                     x_med(1,3) = sign(angle_combs(2,h)-(pi/2)); 
%                 end
%                 
        otherwise
                 x_med(1,1) = cos(theta_med(1));
                 x_med(1,2) = sin(theta_med(1));
                 
                 for c = 2:(n-1)
                     if abs(theta_med(c)-(pi/2)) ~= (pi/2)
                     x_med(1,c+1) = x_med(1,c)*tan(theta_med(c)-(pi/2));                
                     else
                        for cc = 1:c                            
                            x_med(1,cc) = 0;
                        end

                     x_med(1,c+1) = sign(theta_med(c)-(pi/2)); 
                     
                     end
                 end
    end
    
    % Normalize middle point
    %x_med = x_med/norm(x_med);
    
    medpoints{i,1} = x_med;
    medpoints{i+(q/2),1} = -x_med; % Apply symmetry for other half of state-space!!!
    
    x_med = [x_med; zeros(1,n)];
    x_med2 = [-x_med; zeros(1,n)]; % Apply symmetry for other half of state-space!!!
    
    
%     % Normalize vertices
%     for r = 1:size(V,1)
%         V(r,:) = V(r,:)/norm(V(r,:));
%     end
    
    V = [V; zeros(1,n)];    % Add origin as a vertex
    
    Polyhedra{i,1} = Polyhedron('V',V);
    Polyhedra{i+(q/2)} = Polyhedron('V',-V); % Apply symmetry for polyhedra in other half of state-space!!!
    
    % Plot polyhedra
    if n == 3
    plot(Polyhedra{i,1})
    plot(Polyhedra{i+(q/2),1})
        if mod(i,2) == 1 
        plot(Polyhedra{i,1},'color','b')
        plot(Polyhedra{i+(q/2),1},'color','b')
        end
    plot3(x_med(:,1)',x_med(:,2)',x_med(:,3)','k','linewidth',3)
    plot3(x_med2(:,1)',x_med2(:,2)',x_med2(:,3)','k','linewidth',3)
    end
    
 
           
end

%% Check middle point regions
poly_regions = [];  % To store region numbers for each polyhedron ACCORDING TO REG_DET_Q!!!

for i = 1:q
    
poly = Polyhedra{i,1};  % For polyhedron i...                         
center = medpoints{i,1};
%... find region number according to reg_det_Q.
% This index is needed to find the inter-sample time bounds for this region!!!
region_poly = reg_Det_Q(center,m); 
poly_regions = [poly_regions; region_poly];

center_line = [center; zeros(1,n)];

%plot3(center_line(:,1)',center_line(:,2)',center_line(:,3)','k')
%plot(poly)
    
end

U = unique(poly_regions);
unique_regions = size(U,1);
number_regions = q;
% unique should equal number_regions!!!

%% Reachability Analysis MPT

% Because of symmetry, time bounds in the 2nd half of the state space are 
% the same as in the first half. Therefore s-th region has the same
% sampling time bounds as the (s+q/2)-th region
Tau_s_opt = [Tau_s_opt Tau_s_opt];  
Tau_s_max = [Tau_s_max Tau_s_max];
%toc;
INT_POINTS = [];
ReachableSets = cell(q,1);

Reachable_regions = cell(q,1);          % To store reachable region (numbering according to Polyhedra) 
Reachable_regions_regDetQ = cell(q,1);  % To store reachable region (numbering according to reg_Det_Q)
Region_intersections = cell(q,q);
has_region_intersected = zeros(q,q);
%global p
p = 1;
D = parallel.pool.DataQueue;
afterEach(D, @updateWaitBar);
fprintf("finished reachability analysis %3.0f%%\n", 0)
parfor i = 1:q % Loop on region polyhedra

region_reachability = i;
%sprintf(i)

R0 = Polyhedra{i,1};                   % Domain of states (the current region)
timebound_index = poly_regions(i);     % Find region number according to reg_Det_Q!!!
tau_min = Tau_s_opt(timebound_index);  % Start time, which is the region's lower sampling time bound
tau_max = Tau_s_max(timebound_index);  % Final time, which is the region's upper sampling time bound    

h = 10;
Ts = tau_min/h;                        % Define sampling time for reachability analysis
N_min = floor(tau_min/Ts);             % Number of steps for the lower sampling time bound
N_max = ceil(tau_max/Ts);              % Number of steps for the upper sampling time bound
%fprintf('\b\b\b\b\b%3.0f%%\n', 100*(i)/(q)) 
   for step = N_min:N_max
        sysd = c2d(ss(A,B,zeros(1,n),0),step*Ts);                   % Discretize the system
        sys = LTISystem('A',(sysd.A-(sysd.B)*K),'Ts',tau_max);      % Define system for reachability analysis
        R = sys.reachableSet('X',R0,'N',1,'direction','forward');   % Compute reachable set
            
        for j = 1:q     % Loop on regions
            if (has_region_intersected(i,j) == 1)
                continue
            end
            intersection = intersect(R, 100*Polyhedra{j,1});     % Find intersection between R (reachable set) and j-th region polyhedron
            Region_intersections{i,j} = intersection;
                
            if intersection.isEmptySet == 0 && size(intersection.V,1) > 1   % If intersection is not empty and not just a point (only 1 vertex)...

                % Region transitions with region numbering according to
                % polyhedra numbering
                Reachable_regions{i,1} = [Reachable_regions{i,1}, j];       % ... add j-th region to reachable regions of region i
                has_region_intersected(i,j) = 1;
                % Region transitions with region numbering according to
                % reg_det_Q!!!
%                 Reachable_regions_regDetQ{poly_regions(i),1} = [Reachable_regions_regDetQ{poly_regions(i),1}, poly_regions(j)];


            elseif intersection.isEmptySet == 0 && size(intersection.V,1) == 1  &&  norm(intersection.V,1) >= 1e-5 % If intersection is not empty and not just a point (1 vertex) that is the origin...

                % Region transitions with region numbering according to
                % polyhedra numbering
                Reachable_regions{i,1} = [Reachable_regions{i,1}, j];       % ... add j-th region to reachable regions of region i
                has_region_intersected(i,j) = 1;
                % Region transitions with region numbering according to
                % reg_det_Q!!!
%                 Reachable_regions_regDetQ{poly_regions(i),1} = [Reachable_regions_regDetQ{poly_regions(i),1}, poly_regions(j)];
            end

%                     if size(intersection.V,1) == 1
%                         INT_POINTS = [INT_POINTS; intersection.V]; 
%                     end
                
                
                    %end
            end
    end

%fprintf('\b\b\b\b\b%3.0f%%\n', 100*(i)/(q)) 
send(D,i);

end
%% convert polyhydra regions to regDetQ regions
for i = 1:q
    Reachable_regions_regDetQ{poly_regions(i),1} = arrayfun(@(j) poly_regions(j), Reachable_regions{i,1});
end


%% Filter out double regions
for ii = 1:q
    Reachable_regions{ii,1} = unique(Reachable_regions{ii,1});
    Reachable_regions_regDetQ{ii,1} = unique(Reachable_regions_regDetQ{ii,1});
end

function updateWaitBar(~)
    fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*p/q))
    p = p+1;
end

%% CHECK REGION TRANSITIONS USING POLYHEDRA
%  Find region sequence with numbering according to the cell Polyhedra,
%  instead of the numbering according to reg_Det_Q.
%  (Optional, just to double check)

% reg_Seq_Polyhedra = [];
% x_check = [x_0, x_trig];
% 
% for k = 1:size(medpoints,1)   %(x_check,2);
%     
%     %x = x_check(:,k);
%     %
%     x = medpoints{k,1}';
%     x = Polyhedron('V',x');
%     
%     for i = 1:q % Loop on polyhedra
%         state_intersect = intersect(x, 100*Polyhedra{i,1});
%         if state_intersect.isEmptySet == 0
%             reg_Seq_Polyhedra = [reg_Seq_Polyhedra,i];
%             %break
%         end
%     end
%     
% end
% %% Check for impossible transitions (transitions not predicted to be possible by reachability analysis
% 
% impossible_transitions = 0;
% %trouble_poly = [];             
% % To save polyhedra for which reachability analysis failed; requires the
% % previous section to be uncommented!!!
% 
% for i = 1:(size(reg_Seq,2)-1)
%     region = reg_Seq(i);
%     region_next = reg_Seq(i+1);
%     
%     % If a region is not among the reachable regions for the previous
%     % region, the transition was not predicted to be possible by the
%     % reachability analysis
%     if ismember(region_next, Reachable_regions_regDetQ{region,1}) == 0
%         impossible_transitions = impossible_transitions+1;
%         %trouble_poly = [trouble_poly, reg_Seq_Polyhedra(i)];   
%     end
% end
% 
% % Polyhedron for which reachability analysis is not correct
% %trouble_poly = unique(trouble_poly);
% 
% % Number of impossible transitions (should be zero)
% impossible_transitions
% 
% endtime_reach = toc
end

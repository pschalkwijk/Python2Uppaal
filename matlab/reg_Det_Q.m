%% Region determination of a state
% Note: if a state lives on the boundaries of two intersection then the
% output vector y will have 2^(n-1) positive entries.

function x_region = reg_Det_Q(x,m)

n = length(x);

%% Initialize regions for state-space abstraction
Q = Q_2_fun(2,m); % generate Q matrices for two-dimensional plane, given the number of subdivions m

% Obtain nD regions by combining 2D regions (projections of nD regions on (n-1) 2D planes)
Two_dim_div = cell((n-1),1);

for i = 1:(n-2)
    Two_dim_div{i,1} = [1:m];                   % ((n-2) x m - 2D regions
end

Two_dim_div{n-1,1} = [1:2*m];                   % 2*m - 2D regions

AllRegions = combvec(Two_dim_div{1:(n-1),1});   % All possible combinations of 2D regions (each combination corresponding to an nD region)

reg_seq_2D = zeros((n-1),1);          % To save sequence of found regions
reg_seq_2D_corrQ = zeros((n-1),1);    % To save sequence of EXACT CORRECT Q MATRICES ('mirroring' for ALL 2D planes)

%% Check projections onto two-dimensional planes
for i = 1:(n-1)
    
   x_check = [x(i); x(i+1)];                    % Pick two successive indices for projection onto two-dimensional plane
         
   for j=1:m                                    % For all possible regions in a two-dimensional plane...
       if  x_check'*Q{j,1}*x_check >= 0          % ... check whether x_check lies in this region
           region_2D = j;
           region_2D_corrQ = j;
           

            if  i==n-1 && x_check(2) < 0   % Making use of symmetry, mirror region if last dimension is negative
                region_2D = region_2D + m;
            end
            
            if  x_check(2) < 0 % mirroring applied on 2D planes REGARDSLESS OF i !!!
                region_2D_corrQ = region_2D_corrQ + m;
            end                   

                reg_seq_2D(i,1) = region_2D;         % Save found 2D region
                reg_seq_2D_corrQ(i,1) = region_2D_corrQ;
           break
       end
   end
     
end

%% Find region for x

for i = 1:size(AllRegions,2)
    
    if min(reg_seq_2D == AllRegions(:,i)) == 1;     % Check to which Region (column of AllRegions) reg_seq_2D (containing 2D proj. indices for x) corresponds
        x_region = i;
        %break
    end
end
end
    

     
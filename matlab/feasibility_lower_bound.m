function feas = feasibility_lower_bound(tau, l, sigma, N, n, L, Q, region, ops, tol, nu, m, f)
% Region increment
Phi = Phi_fun_etc(tau,l,sigma,N,L);
%         clear con_break
con_break = 0;                        % Loop break condition
if f==1
    for y=0:floor(tau*l/sigma)
        % increment in each row of Phi
        for z=1:N+1                   
            % increment in each column of Phi
            clear U
            e = sdpvar((n-1),1);  % generate (n-1) scalars e(psilon)
            Q_found = cell((n-1),1);
            % Read Q-matrices belonging to current region
            for i = 1:(n-1)
                index_Qfound = region(i,1);              
                Q_found{i,1} = Q{index_Qfound,1};     
            end 
            % Convert 2x2 Q-matrices to nxn matrices by adding zeros at indices of states that are NOT considered
            Q_LMI = cell((n-1),1);

            for i = 1:(n-1)     
                Q_2D = Q_found{i,1};
                Q_nD = zeros(n,n);
                Q_nD(i,i) = Q_2D(1,1);
                Q_nD(i,i+1) = Q_2D(1,2);
                Q_nD(i+1,i) = Q_2D(2,1);
                Q_nD(i+1,i+1) = Q_2D(2,2);
                Q_LMI{i,1} = Q_nD; 
            end                            
            inEq = 0;
            for i = 1:(n-1)                               
                inEq = inEq + e(i)*Q_LMI{i,1};    % Build up LMI from combination of 2D Q-matrices and scalars e(psilon)    
            end
            Con_e = [e(:)>=tol, (Phi{z,y+1}+nu*eye(n)+inEq)<=-10^(-5)];
            diag_sol = solvesdp(Con_e,[],ops); 
            %epsilon_lower = [epsilon_lower, double(e)]; 
            % In case of an unacceptable solution from solver
            if diag_sol.problem ~= 0
                con_break = 1;
                yalmiperror(diag_sol.problem);
                break;         % Acting on z
            end
        end     
        if con_break == 1      % Acting on y
            break;
        end
    end                        % y loop-on (s'th) horizontal blocks
else
    for y=(floor(tau*l/sigma)-1):...
        floor(tau*l/sigma)
        for z=1:N+1                  
            clear U
            e = sdpvar((n-1),1);  % Generate (n-1) scalars e(psilon)
            Q_found = cell((n-1),1);
            %Read Q-matrices belonging to current region
            for i = 1:(n-1)
                index_Qfound = region(i,1);
                if index_Qfound > m
                    index_Qfound = index_Qfound - m;
                end                                
                Q_found{i,1} = Q{index_Qfound,1};    
            end

            % Convert 2x2 Q-matrices to nxn matrices by adding zeros at indices of states that are NOT considered
            Q_LMI = cell((n-1),1);
            for i = 1:(n-1)                         
                Q_2D = Q_found{i,1};
                Q_nD = zeros(n,n);
                Q_nD(i,i) = Q_2D(1,1);
                Q_nD(i,i+1) = Q_2D(1,2);
                Q_nD(i+1,i) = Q_2D(2,1);
                Q_nD(i+1,i+1) = Q_2D(2,2);
                Q_LMI{i,1} = Q_nD;                         
            end

            inEq = 0;
            for i = 1:(n-1)
                inEq = inEq + e(i)*Q_LMI{i,1}; % Build up LMI from combination of 2D Q-matrices and scalars e(psilon)
            end                            
            Con_e = [e(:)>=tol, (Phi{z,y+1}+nu*eye(n)+inEq)<=-10^(-5)];
            diag_sol = solvesdp(Con_e,[],ops);                   
            e_ks = double(e);
            if diag_sol.problem ~= 0
                con_break = 1;
                yalmiperror(diag_sol.problem);
                break;
            end
        end                         
        if con_break == 1            
            break;
        end
    end 
end
feas = con_break;
end
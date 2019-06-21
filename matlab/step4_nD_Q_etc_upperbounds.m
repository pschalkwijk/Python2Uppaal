%% Deriving the regional maximum sampling period. 
% This script calculates lower time bounds on the sampling time for each
% region for an nD system using projections of each region onto 2D planes
% and their corresponding 2x2 Q-matrices.

% Time step parameters preferably already defined in a control system
% specific run script!!!
%sigma_max = 2; % --> this is a design parameter for an estimation of the
                % upper approximation of the regional time period  

% l_star = sigma_max/0.01; %--> Note: 0.01 is related to the desired
%                               precision for convex embedding
% del_tau_4 = 0.001;            %step by step change this one
function Tau_s_max = step4_nD_Q_etc_upperbounds(m, l_star, N_conv, alpha, ...
                                            A, B, K, n, Tau_s_opt, del_tau_4, ...
                                            AllRegions, Q, sigma_max, nu, eps, tol)
    q = m^(n-1); 
    % Note: here l could be changed into l_star for better accuracy
    LL = LL_fun_etc(sigma_max,l_star,N_conv,alpha,A,B,K);  
    
    ops = sdpsettings('solver','sedumi','sedumi.eps',eps,...
        'sedumi.maxiter',100,'verbose',0, 'cachesolvers', 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Forward search from Tau_s_opt to find the first feasible time

    Tau_s_max = Tau_s_opt; 
    % Note: based on the forward feasibility search;
    % Tau_s_opt contains lower bound tau for each region

    p = 1;
    D = parallel.pool.DataQueue;
    afterEach(D, @updateWaitBar);
    fprintf('finished upper regions: %3d%%\n', 0)
    clear mm ff y z primalfeas1 dualfeas1
    
    parfor mm=1:q                                    % Region increment
        
        region = AllRegions(:,mm);
        tau_s_opt = Tau_s_opt(1,mm);
        
        a = tau_s_opt;
        c = sigma_max;
        % check beginning point for gss
        xa = feasibility_upper_bound(a, l_star, sigma_max, ...
            N_conv, n, LL, Q, region, ops, tol, nu);
        xc = feasibility_upper_bound(c, l_star, sigma_max, ...
            N_conv, n, LL, Q, region, ops, tol, nu);
        f = (sigma_max-tau_s_opt)/del_tau_4
        if or(xc ~= 0, xa ==0)
%             disp('sigma max should be a feasible solution')
            % In the original code, tau_s_max = tau_s_opt in this case.
            Tau_s_max(1, mm) = tau_s_opt;
        else
            while(abs(c-a)>del_tau_4)
                b = (c+a)/2;
                xb = feasibility_upper_bound(b,l_star,sigma_max, N_conv,...
                     n, LL, Q, region, ops, tol, nu);
                if xb == 0
                    c = b
                else
                    a = b
                end
            end
            Tau_s_max(1,mm) = a
        end
        send(D, mm);
    end
    figure;
    hold on
    title('Upper bounds on sampling time')
    xlabel('Region')
    ylabel('\tau_{s}^{v}')
    plot(Tau_s_max, '*')
    grid on
    drawnow
    
    function updateWaitBar(~)
        fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*p/q))
        p = p+1;
    end
end





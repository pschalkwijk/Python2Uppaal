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
function Tau_s_max = step4_nD_Q_etc_upperbounds_gss(m, l_star, N_conv, alpha, ...
                                            A, B, K, n, Tau_s_opt, del_tau_4, ...
                                            AllRegions, Q, sigma_max, nu, eps, tol)
    q = m^(n-1); 
    % Note: here l could be changed into l_star for better accuracy
    LL = LL_fun_etc(sigma_max,l_star,N_conv,alpha,A,B,K);  
    
    ops = sdpsettings('solver','sedumi','sedumi.eps',eps,...
        'sedumi.maxiter',100,'verbose',0, 'cachesolvers', 1);
    
    s = (1 + sqrt(5))/2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Forward search from Tau_s_opt to find the first feasible time

    Tau_s_max = Tau_s_opt; % Note: based on the forward feasibility search; Tau_s_opt contains lower bound tau for each region

    epsilon_upper = [];
    p = 1;
    D = parallel.pool.DataQueue;
    afterEach(D, @updateWaitBar);
    fprintf('finished upper regions: %3d%%\n', 0)
    clear mm ff y z primalfeas1 dualfeas1
    parfor mm=1:q                                    % Region increment
        region = AllRegions(:,mm);
        tau_s_opt = Tau_s_opt(1,mm);
        % f is the number of tau_star considered
        f = floor((sigma_max-tau_s_opt)/del_tau_4); 
        a = tau_s_opt;
        d = sigma_max;
        % check beginning point for gss
        xa = feasibility_upper_bound(a, l_star, sigma_max, ...
            N_conv, n, LL, Q, region, ops, tol, nu);
        xd = feasibility_upper_bound(d, l_star, sigma_max, ...
            N_conv, n, LL, Q, region, ops, tol, nu);
        
        if xd ~= 0
            % disp('sigma max should be a feasible solution')
            Tau_s_max(1, mm) = sigma_max;
        end
        
        if (xa ~= 0) && (xd == 0) 
            b = a;
            c = d;
            while (abs(a - d) >= del_tau_4)
                b = a + (d-a)/s;
                c = d - (d-a)/s; 
                if b > c
                    [b, c] = deal(b,c);
                end
                
                fb = feasibility_upper_bound(tau_s_opt+b, l_star, sigma_max,...
                    N_conv, n, LL, Q, region, ops, tol, nu);
                fc = feasibility_upper_bound(tau_s_opt+c, l_star, sigma_max,...
                        N_conv, n, LL, Q, region, ops, tol, nu);
                if (fc == 0)
                    d = c;
                else
                    b = c;
                    fb = fc;
                end
                
                if (fb == 0)
                    d = b;
                else
                    a = b;
                end
            end
            Tau_s_max(1,mm) = tau_s_opt+d;
        end
        
%     fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*mm/q))
        send(D, mm);
    end
    figure
    plot(Tau_s_max,'*', 'Color', 'b'),grid on
    title('Upper bounds on sampling time')
    xlabel('Region')
    ylabel('\tau_{s}^{v}')
% 
    function updateWaitBar(~)
        fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*p/q))
        p = p+1;
    end
end





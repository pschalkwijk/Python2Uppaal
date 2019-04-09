%% Regional lower bound approximation for n-dimension with diagonal approach
% This script calculates lower time bounds on the sampling time for each
% region for an nD system using projections of each region onto 2D planes
% and their corresponding 2x2 Q-matrices.
function Tau_s_opt = step3_nD_Q_etc_lowerbounds_gss(m, sigma_bar, l, N_conv, ...
                                alpha, A, B, K, n, tau_opt_nu, del_tau_3, ...
                                AllRegions, Q, nu, eps, tol)
%%
q = m^(n-1);
% Time step parameters preferably already defined in a control system
% specific run script!!!
% del_tau_3 = sigma_bar/100;           

LL = LL_fun_etc(sigma_bar,l,N_conv,alpha,A,B,K);

ops = sdpsettings('solver','sedumi','sedumi.eps',eps,...
    'sedumi.maxiter',100,'verbose',0, 'cachesolvers', 1);

Tau_s_opt = tau_opt_nu*ones(1,q);              % Sampling function vector

epsilon_lower = [];

% create a progress bar
p = 1;
D = parallel.pool.DataQueue;
afterEach(D, @updateWaitBar);
clear mm ff y z primalfeas1 dualfeas1
f = floor((sigma_bar-tau_opt_nu)/del_tau_3); 
fprintf('finished lower regions: %3.0f%%\n', (100.0*0/q))
parfor mm=1:q   % Region increment 
%     mm = 1;
    % Find combination of Q-matrices corresponding to current region
    Q_sec = AllRegions(:,mm);
    %f is the number of tau_star considered
    a = 0;
    fa = feasibility_lower_bound(tau_opt_nu, l, sigma_bar, N_conv, n, LL,...
        Q, Q_sec, ops, tol, nu, m, f);
    d = sigma_bar;
    fd = feasibility_lower_bound(sigma_bar, l, sigma_bar, N_conv, n, LL,...
        Q, Q_sec, ops, tol, nu, m, f);
    
    s = (1 + sqrt(5))/2;
    if (fa ~= 0)
%         disp('lower bound satifies conditions')
        % Use the lower bound
    elseif (fd == 0)
        disp("upper bound doesn't satisfiy conditions")
        % line search down to find appropriate upper bound
    elseif ((fa == 0) && (fd ~= 0))
%         disp("starting linesearch")
        b = a;
        c = d;
        while (abs(a-d) >= del_tau_3)
%             disp(b-c);
            b = a + (d - a)/s;
            c = d - (d - a)/s;
            if b > c
                [b, c] = deal(c,b);
            end
            fb = feasibility_lower_bound(tau_opt_nu+b, l, sigma_bar, N_conv, n, LL,...
                Q, Q_sec, ops, tol, nu, m, f);
            fc = feasibility_lower_bound(tau_opt_nu+c, l, sigma_bar, N_conv, n, LL,...
                Q, Q_sec, ops, tol, nu, m, f);
            if (fb == 0)
                a = b;
            else
                c = b;
                fc = fb;
            end
            if (fc == 0)
                a = c;
            else
                d = c;
            end
            Tau_s_opt(1,mm) = tau_opt_nu + a;
        end        
    else
        sprintf("bounds not as expected")
    end
    send(D, mm);
%     fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*mm/q))
end   
   
% Check if all epsilon's are positive!!!
% min_epsilon_lower = min(min(epsilon_lower));
% 
% figure 
% plot(Tau_s_opt,'*')


figure
plot(Tau_s_opt,'*', 'Color','b')
grid on
title('Lower bounds on sampling time')
xlabel('Region')
ylabel('\tau_{s}^{v}')

function updateWaitBar(~)
    fprintf('\b\b\b\b\b%3.0f%%\n', (100.0*p/q))
    p = p+1;
end

end
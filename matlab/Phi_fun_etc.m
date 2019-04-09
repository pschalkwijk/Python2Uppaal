function Phi = Phi_fun_etc(tau_s,l,sigma_bar,N_conv,LL)
% Calculating Phi (i,j) = Phi_hat (i,j); see eq.14 in 
% 'A State Dependent Sampling for Linear State Feedback' by Fiter, Hetel,
% Perruquetti and Richard. Equation references in this function refer to 
% this paper.

% i in {0,...,N_conv} & j in {0,...,floor(tau_s*l/sigma_bar)}
% Note: here the Nu is not yet considered in the Phi's structure!

sigma_prim = sigma_bar/l;
if tau_s<sigma_bar
 Phi = cell(N_conv+1,floor(tau_s*l/sigma_bar)+1);  
else
 Phi = cell(N_conv+1,floor(tau_s*l/sigma_bar));
end

if tau_s<sigma_bar
    for j=0:floor(tau_s*l/sigma_bar)     % j in Lemma 6
        if j<floor(tau_s*l/sigma_bar)
            Phi{1,j+1} = LL{1,j+1};
            for w=2:N_conv+1                % w as (i-1) in eq. 15
                Phi{w,j+1} = Phi{w-1,j+1}+sigma_prim^(w-1)*LL{w,j+1};
            end

        else
            Phi{1,j+1} = LL{1,j+1};
            for w=2:N_conv+1                % w as (i-1) in eq. 15
                Phi{w,j+1} = Phi{w-1,j+1}+(tau_s-j*sigma_prim)^(w-1)*LL{w,j+1};
            end
        end

    end
else
    for j=0:(floor(tau_s*l/sigma_bar)-1)    % Note the difference with the above j-counter
        Phi{1,j+1} = LL{1,j+1}; 
        for w=2:N_conv+1
%             sigma_prim^(w-1)% w as (i-1) in eq. 15
            Phi{w,j+1} = Phi{w-1,j+1}+sigma_prim^(w-1)*LL{w,j+1};
        end
    end
end
end









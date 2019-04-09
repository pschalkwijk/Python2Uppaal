%% Considering exponential stability-solving LMI's (x^T) * Phi(i,j) * x <=0  
% This script contains a line search on tau_s to find the minimal value
% for it (over the entire state space).

%% Design parameters
%  Time step parameters preferably already defined in a control system
%  specific run script!!!
%  tau_min = 0;                   
%  del_tau = sigma_bar/100;

gg = floor((sigma_bar-tau_min)/del_tau);  % Number of steps considered for tau_s

tau_opt = 0;
LL = LL_fun_etc(sigma_bar,l,N_conv,alpha,A,B,K);
% size(LL);
% LL{:,1};

for v=0:gg

    clear eig_Phi_max
    tau_s = tau_min+v*del_tau;
    Phi = Phi_fun_etc(tau_s,l,sigma_bar,N_conv,LL);
%     size(Phi)
    eig_Phi_max = zeros(size(Phi));
    for j=0:floor(tau_s*l/sigma_bar)
            for i=1:N_conv+1
                eig_Phi_max(i,j+1) = max(double(max(eigs(Phi{i,j+1}))));
            end
    end
    
    if max(max(eig_Phi_max))<=0
        tau_opt = tau_s;
    else
        break;
    end
     
end

%% Checking whether the resulting P and tau_min are also valid numerically

% LL = LL_fun_etc(sigma_bar,l,N_conv,alpha,A,B,K);
% Phi = Phi_fun_etc(tau_opt,l,sigma_bar,N_conv,LL);
% 
% eig_Phi_max = zeros(size(Phi));
% for j=0:floor(tau_opt*l/sigma_bar)
%         for i=1:N_conv+1
%             eig_Phi_max(i,j+1) = max(double(max(eig(Phi{i,j+1}))));
%         end
% end
% 
% % eig_Phi_Max = max(eig_Phi_max); % Should be nonpositive
% for rr=1:N_conv+1
%     figure
%     grid on
%     plot((0:j),eig_Phi_max(rr,:),'*')   % Should all be nonpositive, because Phi should be seminegative definite
%     xlabel('$j$','interpreter','latex')
%     ylabel('$\lambda_{max}\Phi_{(1:N+1,j)}$','interpreter','latex')
% end
%     
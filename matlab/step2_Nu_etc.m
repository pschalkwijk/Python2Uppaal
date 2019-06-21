%% Nu determination 
% This script calculates the value for Nu as described in the paper
% 'A State Dependent Sampling for Linear State Feedback' by Fiter, Hetel,
% Perruquetti and Richard. Equation references in this script refer to this
% paper.

% del_sig = 0.001;             % sigma_prime increment --> eq. 19 
% Preferably already defined in a control system
% specific run script!
LL = LL_fun_etc(sigma_bar,l,N_conv,alpha,A,B,K); % eq. 16

yy = floor((sigma_bar/l)/del_sig)+1 ;  % sigma_prime lies in [0,sigma_bar/l]
eig_max = zeros(yy,l);                 % Size form := (sigma_prim,l)
      
% Numerical integration method for matrix exponential calculation
% Able to handle singular matrices
eAt = @(t) expm(A*t);               
%     
int_eAx = @(x) integral(eAt,0,x,'ArrayValued',true);
ABK = A-B*K;
Phi_con = @(x)(eye(n)-(eye(n)+int_eAx(x)*ABK))'*...
                (eye(n)-(eye(n)+int_eAx(x)*ABK))-...
                alpha*(eye(n)+int_eAx(x)*ABK)'*...
                (eye(n)+int_eAx(x)*ABK);      
     

for y=1:yy  % sigma_prime counter
    
    sig_prim = del_sig*(y-1);
    Phi_hat = Phi_hat_fun(l,sig_prim,N_conv,LL);  
    for r=0:l-1                       
        phi_con = Phi_con(sig_prim+r*sigma_bar/l);
        sub_Val = phi_con-double(Phi_hat{1,r+1});  
        ev = eigs(sub_Val);
        eig_max(y,r+1) = max(double(abs(ev)));   
    end
    
end
nu = max(max(abs(eig_max)));
disp(nu);

%% 2.2: Checking that Phi_k <=0 (Nu is added here)
tau_s = tau_opt;                         

% del_tau_2 = 0.001;                         % tau_s increment length
% Preferably already defined in a control system
% specific run script!

clear i j v tau_op_nu con_Break
con_Break = 0;                           % For loop stop condition

for v=tau_opt:-del_tau_2:0               
    tau_val = v;
    Phi = Phi_fun_etc(tau_val,l,sigma_bar,N_conv,LL);
    clear eig_Phi_max_nu
    eig_Phi_max_nu = zeros(size(Phi));
    for j=0:floor(tau_val*l/sigma_bar)
        for i=1:N_conv+1
             eig_Phi_max_nu(i,j+1) = max(double(eigs(Phi{i,j+1}+nu*eye(n))));
        end
    end
    if max(max(eig_Phi_max_nu))<=0
        tau_opt_nu = tau_val;
        break;                           
    end
end

%% Plotting the maximal eigenvalues of Phi

% Phi = Phi_fun_etc(tau_opt_nu,l,sigma_bar,N_conv,LL);
eig_Phi_max_nu = zeros(size(Phi));

for j=0:floor(tau_opt_nu*l/sigma_bar)
    for i=1:N_conv+1
         eig_Phi_max_nu(i,j+1) = max(double(eig(Phi{i,j+1}+nu*eye(n))));
    end
end


figure
hold on
grid on
for r_1=1:N_conv+1
    plot(eig_Phi_max_nu(r_1,:)) 
    xlabel('$j$','interpreter','latex')
    ylabel('$\lambda_{max}\Phi_{(1:N+1,j)}$','interpreter','latex')
end

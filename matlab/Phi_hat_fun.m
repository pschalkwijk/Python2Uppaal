function Phi_hat = Phi_hat_fun(l,sig_prim,N_conv,LL)
%%%% calculating Phi_bar (N,j); eq. 23 in Fiter's paper
% sigma_prim = sigma_bar/l;
Phi_hat = LL(1,:);                     % Phi_(0,j)

for j=0:l-1                            % (j) See Lemma 6 in Fiter's paper
    
    for w=2:N_conv+1                % (w-1=k) in eq. 23
        Phi_hat{1,j+1} = Phi_hat{1,j+1}+sig_prim^(w-1)*LL{w,j+1};
    end

end

end
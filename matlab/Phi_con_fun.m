function Phicon = Phi_con_fun(A,x) % Calculates Phi from eq. 19 in Fiter's paper
        
        n = size(A); % System order
        
% Numerical integration method for matrix exponential calculation
% Able to handle singular matrices
        eAt = @(t) expm(A*t);               
        inteAx = integral(eAt,0,x,'ArrayValued',true);
        
        % calculating Phi_con 
        Phicon =(eye(n)-(eye(n)+inteAx(x)*(A-B*K)))'*...
               (eye(n)-(eye(n)+inteAx(x)*(A-B*K)))-...
               alpha*(eye(n)+inteAx(x)*(A-B*K))'*...
               (eye(n)+inteAx(x)*(A-B*K));      
end
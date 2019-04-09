function MM = M_fun(sigma_prime, A) %calculates Mj (eq.18)
    % Numerical integration method for matrix exponential calculation
    % Able to handle singular matrices
    % A = [0, 1; -2, 3]
    eAt = @(t) expm(A*t);               
    MM = integral(eAt,0,sigma_prime,'ArrayValued',true);
end
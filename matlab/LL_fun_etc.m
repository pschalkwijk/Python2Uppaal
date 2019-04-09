function LL = LL_fun_etc(sigma_bar,l,N_conv,alpha,A,B,K)
% Calculating L-matrices; see eq. 16 in
% 'A State Dependent Sampling for Linear State Feedback' by Fiter, Hetel,
% Perruquetti and Richard. Equation references in this function refer to 
% this paper.
% NOTE: See: A Formal Traffic Characterization Of LTI Event-triggered
% systems for the actual formula's: (13, 14, 15, 17)
% Note: N_conv must be greater than 2
n = length(A(1,:));

    function MM = M_fun(sigma_prime) %calculates Mj (eq.18 or eq.15 (Kolajirani))
        % Numerical integration method for matrix exponential calculation
        % Able to handle singular matrices
        eAt = @(t) expm(A*t);               
        MM = integral(eAt,0,sigma_prime,'ArrayValued',true);
    end

% For sum in LL_(k>=2):sum_in_L is not affected by the change in j (eq.16 or eq.13) 
sum_in_L = cell(N_conv-1,1);                            % For k=2,...,N_conv
sum_in_L{1,1} = eye(n);                                 % For k=2
if N_conv>2
    for k=3:N_conv                                      % For k>=3
        sum_in_L{k-1,1} = A^(k-2)/factorial(k-1);       % For i=1
        for t=2:k-1                                     % Loop on i to k-1
            sum_in_L{k-1,1} = sum_in_L{k-1,1}+(A^(t-1))'/factorial(t)*...
                A^(k-t-1)/factorial(k-t);                         
        end
    end
end
%%%
LL = cell(N_conv+1,l);

for j = 0:l-1                % Increment of each time subdivision-Note: only up to floor(tau_s*l/sigma_bar) matters based on Lemma 6.
    
    clear sigma_prime M PI_1 N PI_2 ii
    sigma_prime = j*sigma_bar/l;           % Needed to calculate Mj
    M = M_fun(sigma_prime);                % M_j
    PI_1 = eye(n)+M*(A-B*K);               % PI_(1,j)
    N = A*M+eye(n);                        % N_j
    PI_2 = N*(A-B*K);                      % PI_(2,j)
%     abk = (A-B*K)
    % Calculating [L_(0,j);L_(1,j)]
    LL{1,j+1} = eye(n)-PI_1-PI_1'+(1-alpha)*(PI_1')*PI_1;           %L_{0,j}
    LL{2,j+1} = ((1-alpha)*PI_1'-eye(n))*PI_2+(((1-alpha)*PI_1'-eye(n))*PI_2)';  %L_{1,j}
    
    % Calculating [L_(2,j);...;L_(N_conv+1,j)]
    for ii=2:N_conv                                    % Loop on k>=2 (eq. 16 or 13)
        LL{ii+1,j+1} = ((1-alpha)*PI_1'-eye(n))*(A^(ii-1))/factorial(ii)*...
            PI_2+(((1-alpha)*PI_1'-eye(n))*(A^(ii-1))/factorial(ii)*...
            PI_2)'+(1-alpha)*PI_2'*sum_in_L{ii-1,1}*PI_2;
    end
                  
    
end

end



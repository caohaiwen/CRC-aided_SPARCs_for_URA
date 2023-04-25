function x = FWHT_user_subsampling(mode, beta, n_2power, M, ordering)
% compute a general x = A*beta via a fast walsh-hadamard transform.

% ordering : the rows randomly chosen from Hadamard matrix.

switch (mode)
    
    % mode = 1 , for A*beta 
    case 1 
        beta_2power = zeros(n_2power, 1);
        beta_2power(n_2power - M + 1 : n_2power) = beta;
        xx = fwht_user(beta_2power')';
        x = xx(ordering);
        
    % mode = 2, for A'*z
    case 2
        beta_2power = zeros(n_2power, 1);
        beta_2power(ordering) = beta;
        xx = fwht_user(beta_2power')';
        x = xx(n_2power - M + 1 : n_2power);
end

end


function [PUPE, average_num_checks] = UMAC_lossless_stitching_tree_codes(num_active, num_round, M, r,  num_iterations)
%  function [p_md,p_fa] = unsourcedSPARC(L, S, J, K_a, M, data_profile, P, iter, fading)
%unsourcedSPARC Compute the error rate of the MU-SPARC coding scheme
%
%   num_round:              Number of slots
%   num_active:            active users
%   M:              num of bits per section
%   r:   vector of the number of parity bits per section
%   num_iterations:           number of simulation runs
%
%   struct fading with the following fields
%   'type' = {'no_fading'(default),'uniform','exp','pathloss','shadowing_pathloss'}
%   'lower_limit'(default = 10),'upper_limit'(default = 30): limits for the uniform case

J = log2(M);
N               = 2^J*num_round;
data_profile  = J - r;


% Lookup table for outer code
LT  = createLT(2^max(data_profile));

num_checks = zeros(1, num_iterations);
p_md = 0;
p_fa = 0;

for iter=1:num_iterations
    
    % random messages
    msg = generateMessages(data_profile, num_active);

    % C is the matrix of parity checks needed for the decoder
    [symbol, C] = outerEncoder( msg, 2^J, num_round, data_profile );

    for num_r = 1 : num_round
        symbol_list{num_r} = symbol(num_r, :)-1;
    end
    
    [msg_list, num_checks(iter)]    = outerDecoder(symbol_list, C, r, 2^J, LT);

    [fn,fp]     = countErrors(msg_list, msg);
    p_md        = p_md + fn/num_active;
    p_fa        = p_fa + fp/length(msg_list);
end
p_md = p_md/num_iterations;
p_fa = p_fa/num_iterations;
PUPE = p_md + p_fa;
average_num_checks = sum(num_checks) / num_iterations;

end



%create lookuptable for sums of bits
function LT = createLT(N)
    LT = zeros(N,1);
    for i = 1:N
        LT(i) = mod(sum(dec2bin(i-1)-'0'),2);
    end
end

% Return: LxK matrix of symbols where each symbol has data_profile(l) bits
function msg_matrix = generateMessages(data_profile, K)
    L           = length(data_profile);
    msg_matrix  = zeros(L, K);
    for i = 1:L
        msg_matrix(i,:) = randi([0,2^data_profile(i) - 1], 1, K);
    end
end


% Create list of position indices from a support vector
function symbols = getSymbolList(rpos_matrix)
    L = size(rpos_matrix,2);
    for i = 1:L
        % outer decoder wants symbols from 0:B-1;
        symbols{i} = find(rpos_matrix(:,i)==1)' - 1;
    end
end

% Count the messages which miss in the output list, and the messages
% which were not transmitted but appeared in the output list
function [errs, false_positives] = countErrors(est, sparc_matrix)
    K       = size(est,2);
    errs    = size(sparc_matrix,2);
    
    false_positives = 0;
    for i = 1:K
        codeword = est(:,i);
        [~,indx]=ismember(codeword',sparc_matrix','rows');
        if (indx~=0)
            sparc_matrix(:,indx) = [];
            errs = errs -1;
        else
            false_positives = false_positives + 1;
        end        
    end
end


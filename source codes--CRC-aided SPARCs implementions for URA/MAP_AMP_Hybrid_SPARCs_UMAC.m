function decoded_candidates = MAP_AMP_Hybrid_SPARCs_UMAC(y, n, M, K, extra_candidates,ordering, P, T)
%  this function corresponds to our inner decoder which consists of a
%  successive cancellation algorithm and a simplified AMP decoder.

% y : the received signal
% n : the number of channel uses, i.e., the length of transmitted codeword
% M: the section size of SPARCs
% K: the number of active users
% extra_candidates: the number of extra candidates
% ordering: chosen columns of the Hadamard matrix
% P: average power
% T: the maximum iterations for AMP decoding

allocated_P = sqrt(n*P);

decoded_candidates = zeros(1, K+extra_candidates);

n_2power = M;

% to avoid collisions among different active users, we use
% successive-cancellation-based MAP decoder as a preprocessing procedure.
K_sc = min(round(K/5), 10);
for l = 1 : K_sc
    beta_ = zeros(M, 1);
    inner_prod = FWHT_user_subsampling(2, y, n_2power, M, ordering);
    [~, decoded_info_pos_] = max(inner_prod);
    decoded_candidates(l) = decoded_info_pos_;
    beta_(decoded_info_pos_) = 1;
    y = y - sqrt(P)*FWHT_user_subsampling(1, beta_, n_2power, M, ordering);
end

     
% a simplified AMP decoder

K_remain = K - K_sc;
q = K_remain/M; %sparisity


z = y;
last_tao = 0;
beta = zeros(M, 1);

for t = 1 : T
    tao = sqrt(sum(z.^2)/n);
    if (abs(tao-last_tao) < 1e-6)
        break;
    end
    last_tao = tao;
                                        
    r = allocated_P*beta + FWHT_user_subsampling(2, z, n_2power, M, ordering)/sqrt(n);
     
    denoiser_1_part = exp(-(r-allocated_P).^2 / (2*tao^2)); 
    denoiser_0_part = exp(-r.^2/(2*tao^2));
    
    beta = q*denoiser_1_part ./ ((1-q)*denoiser_0_part+q*denoiser_1_part);
    z = y - sqrt(P)*FWHT_user_subsampling(1, beta, n_2power, M, ordering) ... 
                                       + (z/(tao.^2))*P*sum(beta.*(1-beta));
end

[~, topk_pos] = maxk(beta, K_remain+extra_candidates);
decoded_candidates(K_sc+1:K+extra_candidates) = topk_pos;

end




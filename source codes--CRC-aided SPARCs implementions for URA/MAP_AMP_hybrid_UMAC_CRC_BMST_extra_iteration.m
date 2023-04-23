function [PUPE_MD, PUPE_FA, PUPE_MD_again, PUPE_FA_again] = MAP_AMP_hybrid_UMAC_CRC_BMST_extra_iteration(num_active, num_round, num_info_bits, M, n, P, r, protect_sections, num_trials, memory)
  % encoding and decoding of the proposed concatenated coding system for
  % BMST-CRC codes and SPARCs with an extra iteration for fair comparison.

  % PUPE_MD: PUPE for miss detection;
  % PUPE_FA: PUPE for false alarm;
  % PUPE_MD_again: PUPE for miss detection after one extra iteration;
  % PUPE_FA_again: PUPE for false alarm after one extra iteration;

  % num_active: the number of active users;
  % num_round:  the number of chunks;
  % num_info_bits: a vector consisted with entries showing the number of  info bits for each chunk;
  % M: the size of codebook;
  % n: the length of codewords;
  % P: the averag power constraint;
  % r: a vector consisted with entries showing the number of redundant bits for each chunk;
  % protect_sections: indicate the CRC codeword protecting how many previously consecutive chunks;
  % num_trials: the number of simulation runs;
  % memory: the memory parameter used in BMST-CRC encoding;


poly{4} = [1 0 0 1 1]; % up to length 11

% poly{5} = [1 0 0 1 0 1]; % up to length 26
poly{5} = [1 0 1 0 1 1]; % up to length 10

% poly{6} = [1 0 0 0 0 1 1]; % up to length 57
poly{6} = [1 0 0 0 1 1 1]; % up to (information) length 25

poly{7} = [1 0 1 1 0 1 1 1]; % up to length 56

% poly{8} = [1 0 0 1 0 1 1 1 1]; % up to length 119
% poly{8} = [1 0 0 0 0 0 1 1 1]; % up to length 119
poly{8} = [1 1 1 0 1 0 1 1 1]; % up to length 9

% poly{9} = [1 0 1 0 0 1 0 1 1 1]; % up to length 246
% poly{9} = [1 0 1 1 1 1 1 0 1 1]; % up to length 246
% poly{9} = [1 1 0 0 0 0 1 0 1 1]; % up to length 13
poly{9} = [1 0 0 1 1 1 1 0 0 1]; % up to length 8

% poly{10} = [1 1 0 0 0 1 1 0 0 1 1]; % up to length 501
% poly{10} = [1 0 1 0 1 1 1 0 0 1 1]; % up to length 21
% poly{10} = [1 0 1 0 0 0 1 1 1 0 1]; % up to length 12
% poly{10} = [1 0 1 0 0 1 1 0 1 1 1]; % up to length 5

% poly{11} = [1 0 0 1 1 1 1 0 1 0 1 1]; % up to length 4
poly{11} = [1 0 1 0 1 1 1 0 0 0 1 1]; % up to length 12

poly{13} = [1 0 0 0 0 1 0 1 1 0 1 1 1 1]; % up to length 11

poly{14} = [1 0 0 0 1 1 0 1 1 1 0 0 0 1 1]; % up to length 11

L = num_active;
K = num_active;
m = log2(M);

P_l = P;

% randomly generate interleavers
orderings = zeros(memory, m);
inverse_orderings = zeros(memory, m);
for mmr = 1 : memory
  perm_current = randperm(m);
  inverse_perm_current(perm_current) = 1 : m;
  orderings(mmr, :) = perm_current;
  inverse_orderings(mmr, :) = inverse_perm_current;
end


% generate the measurement matrix A by randomly choosing n rows from Hadamard matrix.
n_2power = M;
ordering =  randperm(n_2power-1, n)+1; 

num_info_bits_round(1) = m;

for current_round = 2 : num_round
  num_info_bits_round(current_round) = num_info_bits_round(current_round-1) + m - r(current_round);
end

extra_candidates = 10;
num_candidates = K + extra_candidates;
num_missed_detection_again = zeros(num_trials, 1);
num_false_alarm_again = zeros(num_trials, 1);
num_missed_detection = zeros(num_trials, 1);
num_false_alarm = zeros(num_trials, 1);


for trial = 1 : num_trials
   rng('shuffle');
   outer_encoded = zeros(L, m*num_round);
   message = zeros(L, num_info_bits);

   % outer encoding via CRC codes at each round
   for l = 1 : L
     message(l, :) = randi(2, 1, num_info_bits)-1;
     if (sum(message(l, 1:m))==0) % to avoid taking the first all-one codeword
       message(l, m) = 1;
     end
     outer_encoded(l, :) = outer_encoding_BMST(message(l, :), memory, m, num_round, r, protect_sections, num_info_bits_round, orderings);
   end

   % inner encoding and decoding via SPARCs
   CRC_encoded_pos = zeros(num_round, L);

   decoded_CRC_encoded_candidates = zeros(num_round, num_candidates);

  
   y = zeros(n, num_round);
   
   for current_round = 1 : num_round
     beta = zeros(M, 1);
     last_end = (current_round-1)*m;
     
      % inner encoding via SPARCs patch by patch
     for l = 1 : L
       current_CRC_encoded_bits = outer_encoded(l, last_end+1:last_end+m);
       CRC_encoded_pos(current_round, l) = bi2de(current_CRC_encoded_bits)+1;
       beta(CRC_encoded_pos(current_round, l)) = beta(CRC_encoded_pos(current_round, l)) + 1;
     end
     
     
     x = sqrt(P)*FWHT_user_subsampling(1, beta, n_2power, M, ordering);
    
     z = randn(n, 1);  
     y(:, current_round) = x + z;

     % MAP+AMP hybrid (inner) decoding pitch by pitch
     [decoded_CRC_encoded_candidates(current_round, :)] = MAP_AMP_Hybrid_SPARCs_UMAC(y(:, current_round), n, M, K, extra_candidates, ordering, P_l, 100);

   end
   
      
   % count duplicate roots
   decoded_roots_msg= decoded_CRC_encoded_candidates(1, :);

   num_decoded_duplicate_roots = 0;
   duplicate_decoded_roots = {};
   [~, unique_roots_idx] = unique(decoded_roots_msg, 'stable');
   roots_duplicate_idx = setdiff(1:num_candidates, unique_roots_idx);
   
   while (~isempty(roots_duplicate_idx))
       num_decoded_duplicate_roots = num_decoded_duplicate_roots + 1;
       current_root_idx = roots_duplicate_idx(1);
       current_root_msg = decoded_roots_msg(current_root_idx);
       duplicate_decoded_roots_current = find(decoded_roots_msg==current_root_msg);
       duplicate_decoded_roots{num_decoded_duplicate_roots} = duplicate_decoded_roots_current;
       roots_duplicate_idx = setdiff(roots_duplicate_idx, duplicate_decoded_roots_current);
   end   
   
   % stitching procedure
   final_decoded_users = Outer_UMAC_stitching_CRC_BMST(num_candidates, orderings, decoded_CRC_encoded_candidates, num_round, M, r, protect_sections, memory, poly, duplicate_decoded_roots);

   
   num_decoded_users = size(final_decoded_users, 2);
   
   for c_u = 1 : num_decoded_users
       decoded_msg = final_decoded_users{c_u};
       decoded_roots(c_u) = decoded_msg(1);
   end

   
   % remove possible duplicates
   [~, unique_idx] = unique(decoded_roots, 'stable');
   root_duplicate_idx = setdiff(1:num_decoded_users, unique_idx);
   duplicate_set = [];
   for idx = 1 : length(root_duplicate_idx)
       current = root_duplicate_idx(idx);
       current_user = final_decoded_users{current};
       root_duplicates = setdiff(find(decoded_roots==decoded_roots(current)), current);
       for compare_roots = root_duplicates
           if (~ismember(compare_roots, duplicate_set))
               compare_user = final_decoded_users{compare_roots};
               if (prod(compare_user==current_user))
                   duplicate_set = [duplicate_set, compare_roots];
               end
           end
       end       
   end
   
   if (~isempty(duplicate_set))
       final_decoded_users(duplicate_set) = []; % remove the duplicated copy
       num_decoded_users = size(final_decoded_users, 2);
   end
   
   num_correct_decoded = 0;
   user_roots = CRC_encoded_pos(1, :);
   decoded_msgs = zeros(num_decoded_users, num_round);
   for c_u = 1 : num_decoded_users
       decoded_msg = final_decoded_users{c_u};
       decoded_msgs(c_u, :) = decoded_msg;
       decoded_user = find(user_roots == decoded_msg(1));
       if (decoded_user)
           if (sum(prod(CRC_encoded_pos(:, decoded_user)==decoded_msg'))) % paths with the same roots are included
               num_correct_decoded = num_correct_decoded + 1;
           end
       end
   end
   
   num_missed_detection(trial) = num_active - num_correct_decoded;
   num_false_alarm(trial) = num_decoded_users - num_correct_decoded;
   

   % one extra iteration, i.e., subtract the decoded users' messages and
   % run our proposed decoding algorithm for the remaining;

   num_left = K-num_decoded_users;
   if (num_left ==0)
       num_missed_detection_again(trial)  = num_missed_detection(trial);
       num_false_alarm_again(trial) = num_false_alarm(trial);
       continue;
   end
   
   num_left_candidates = num_left + extra_candidates;
   decoded_CRC_encoded_candidates_again = zeros(num_round, num_left_candidates);
   decoded_CRC_encoded_candidates_again_MRC = zeros(num_round, num_left_candidates);
   for current_round = 1 : num_round
       y_current = y(:, current_round);
       decoded_pos_current = decoded_msgs(:, current_round);
       beta = zeros(M, 1);
       for pos_idx = 1 : num_decoded_users
           beta(decoded_pos_current(pos_idx)) = beta(decoded_pos_current(pos_idx)) + 1;
       end
       x = sqrt(P)*FWHT_user_subsampling(1, beta, n_2power, M, ordering);
       y_current = y_current - x;
       [decoded_CRC_encoded_candidates_again(current_round, :)] = AMP_SPARCs_UMAC(y_current, n, M, num_left, extra_candidates, ordering, P_l, 100);
       [~, decoded_CRC_encoded_candidates_again_MRC(current_round, :)] = maxk(FWHT_user_subsampling(2, y_current, n_2power, M, ordering), num_left_candidates);
   end
   
   % count duplicate roots
   decoded_roots_msg= decoded_CRC_encoded_candidates_again(1, :);

   num_decoded_duplicate_roots = 0;
   duplicate_decoded_roots = {};
   [~, unique_roots_idx] = unique(decoded_roots_msg, 'stable');
   roots_duplicate_idx = setdiff(1:num_left_candidates, unique_roots_idx);
   
   while (~isempty(roots_duplicate_idx))
       num_decoded_duplicate_roots = num_decoded_duplicate_roots + 1;
       current_root_idx = roots_duplicate_idx(1);
       current_root_msg = decoded_roots_msg(current_root_idx);
       duplicate_decoded_roots_current = find(decoded_roots_msg==current_root_msg);
       duplicate_decoded_roots{num_decoded_duplicate_roots} = duplicate_decoded_roots_current;
       roots_duplicate_idx = setdiff(roots_duplicate_idx, duplicate_decoded_roots_current);
   end   
   
   
   % stitching procedure again
   final_decoded_users_again = Outer_UMAC_stitching_CRC_BMST(num_left_candidates, orderings, decoded_CRC_encoded_candidates_again, num_round, M, r, protect_sections, memory, poly, duplicate_decoded_roots);
     
   num_decoded_users_again = size(final_decoded_users_again, 2);
   
   if (num_decoded_users_again == 0)
       num_missed_detection_again(trial)  = num_missed_detection(trial);
       num_false_alarm_again(trial) = num_false_alarm(trial);
       continue;
   end 
   
   
   decoded_roots = [];
   for c_u = 1 : num_decoded_users_again
       decoded_msg = final_decoded_users_again{c_u};
       decoded_roots(c_u) = decoded_msg(1);
   end

   user_roots = CRC_encoded_pos(1, :);
   
   % remove possible duplicates
   [~, unique_idx] = unique(decoded_roots, 'stable');
   root_duplicate_idx = setdiff(1:num_decoded_users_again, unique_idx);
   duplicate_set = [];
   for idx = 1 : length(root_duplicate_idx)
       current = root_duplicate_idx(idx);
       current_user = final_decoded_users_again{current};
       root_duplicates = setdiff(find(decoded_roots==decoded_roots(current)), current);
       for compare_roots = root_duplicates
           if (~ismember(compare_roots, duplicate_set))
               compare_user = final_decoded_users_again{compare_roots};
               if (prod(compare_user==current_user))
                   duplicate_set = [duplicate_set, compare_roots];
               end
           end
       end       
   end
   
   if (~isempty(duplicate_set))
       final_decoded_users_again(duplicate_set) = [];
       num_decoded_users_again = size(final_decoded_users_again, 2);
   end

   
   num_decoded_users_total = num_decoded_users + num_decoded_users_again;

   for c_u_again = 1 : num_decoded_users_again
       decoded_msg = final_decoded_users_again{c_u_again};
       decoded_user = find(user_roots == decoded_msg(1));
       
       if (decoded_user)
           if (sum(prod(CRC_encoded_pos(:, decoded_user)==decoded_msg'))) % paths with the same roots are included
               num_correct_decoded = num_correct_decoded + 1;
           end
       end
       
   end

   num_missed_detection_again(trial) = num_active - num_correct_decoded;
   num_false_alarm_again(trial) = num_decoded_users_total - num_correct_decoded;

    
end

PUPE_MD = sum(num_missed_detection) / (num_trials*num_active);
PUPE_FA = sum(num_false_alarm) / (num_trials*num_active);

PUPE_MD_again = sum(num_missed_detection_again) / (num_trials*num_active);
PUPE_FA_again = sum(num_false_alarm_again) / (num_trials*num_active);

end


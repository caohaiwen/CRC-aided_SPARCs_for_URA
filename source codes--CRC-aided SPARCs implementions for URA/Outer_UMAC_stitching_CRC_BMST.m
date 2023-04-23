function [final_decoded_users] = Outer_UMAC_stitching_CRC_BMST(num_candidates,orderings, outer_decoded_pos, num_round, M, r, protect_sections, memory, poly, duplicate_decoded_roots)
  % tree decoder for UMAC, i.e., stitching procedure
  % only output unordered list of decoded messages transmitted over unsourced MAC.

  % final_decoded_users: the unordered list of decoded messages;

  % num_active: the number of active users;
  % num_round:  the number of chunks;
  % M: the size of codebook;
  % n: the length of codewords;
  % r: a vector consisted with entries showing the number of redundant bits for each chunk;
  % protect_sections: indicate the CRC codeword protecting how many previously consecutive chunks;
  % memory: the memory parameter used in BMST-CRC encoding;
  % poly: CRC generator polynomial;
  % duplicate_decoded_roots: duplicated decoded roots for different decoding trees;

  
  m = log2(M);
  final_decoded_users = [];

   % stitching procedure

   decoded_CRC_encoded_pos = outer_decoded_pos;

   % with the aid of CRC codes, record the survival paths
   
   num_survival_paths = zeros(num_candidates, num_round);
   survival_current_blocks = cell(num_candidates, num_candidates);

   survival_paths = cell(num_candidates, num_round);  % each cell represents a matrix containing all survival paths.
   
   for current_user = 1 : num_candidates
     survival_paths{current_user, 1} = decoded_CRC_encoded_pos(1, current_user);
     num_survival_paths(current_user, 1) = 1;
     survival_current_blocks{current_user, 1} = de2bi(decoded_CRC_encoded_pos(1, current_user)-1, m);
   end


   for current_round = 2 : num_round
       survival_blocks_prev = cell(num_candidates, num_candidates);
       survival_blocks_prev = survival_current_blocks;
       survival_current_blocks = cell(num_candidates, num_candidates);
       for current_root = 1 : num_candidates
           num_paths_current = 0;
           num_paths_prev = num_survival_paths(current_root, current_round-1);
           survival_paths_prev = survival_paths{current_root, current_round-1};
           survival_path_current = [];
           unique_msg_current = unique(decoded_CRC_encoded_pos(current_round, :), "stable");
           num_protect_section = protect_sections(current_round);
           for current_path_num = 1 : num_paths_prev
               survival_path_prev = survival_paths_prev(current_path_num, :);
               previous_blocks = survival_blocks_prev{current_root, current_path_num};
               superimposed_sum = zeros(1, m);
               num_m = 0;
               start_block_BMST_ = max(current_round-memory, 1);
               for prev_batch = current_round-1 : -1 : start_block_BMST_
                   num_m = num_m + 1;
                   superimposed_sum = superimposed_sum + previous_blocks(prev_batch, orderings(num_m,:));
               end
               
               previous_superimposed_sum = mod(superimposed_sum, 2);
               for current_pos_considered = unique_msg_current
                   CRC_poly = poly{r(current_round)};
                   current_path = [survival_path_prev, current_pos_considered];
                   current_blocks = [];
                   [current_blocks,check] = outer_CRC_checks_BMST_efficient(previous_blocks, previous_superimposed_sum, current_path, current_round, num_protect_section, m, memory, r, CRC_poly);
                   if (check)
                       num_paths_current = num_paths_current + 1;
                       survival_path_current(num_paths_current, 1:current_round) = [survival_path_prev, current_pos_considered];
                       survival_current_blocks{current_root, num_paths_current} = current_blocks;
                   end
               end
           end
           num_survival_paths(current_root, current_round) = num_paths_current;
           survival_paths{current_root, current_round} = survival_path_current;
       end
       %         disp(current_round);
   end
   
   
   % peeling decoding procedure in the sense that peel off nodes included in the unique path
   remaining_message = cell(1, num_round);
   for current_round = 1 : num_round
     remaining_message{current_round} = decoded_CRC_encoded_pos(current_round, :);
   end

   num_remainings_last = num_candidates;
   remaining_roots_elimination = 1 : num_candidates;
   final_num_survival_paths = num_survival_paths(:, num_round);
   num_detected_error_trees = sum(final_num_survival_paths==0);
   remaining_roots_elimination(final_num_survival_paths==0) = []; % excluding the detected failure trees.
   num_decoded_users = 0;

while (1)     
        
   while (1)
        num_remainings = 0;
        remaining_roots = [];
        current_remaining_roots_elimination = remaining_roots_elimination;

        % eliminating the corrected messages
        for current_root  = remaining_roots_elimination
          if (num_survival_paths(current_root, num_round)==1)
            flag_error = 0;
            current_pos = zeros(1, num_round);
            current_path = survival_paths{current_root, num_round};
            current_remaining_roots_elimination(current_remaining_roots_elimination==current_root) = [];
            for current_round = 1 : num_round                
                remaining_message_current = remaining_message{current_round};
                corresp_pos = find(remaining_message_current==current_path(current_round));
                size_pos = length(corresp_pos);
                if (size_pos==0)
                    flag_error = 1;
                    break;

                end
                current_pos(current_round) = corresp_pos(1);
            end
            if (~flag_error) % if success, record it and remove it from the remaining set
                num_decoded_users = num_decoded_users + 1;
                final_decoded_users{num_decoded_users} = current_path;
                for current_round = 1 : num_round
                    remaining_message_current = remaining_message{current_round};
                    remaining_message_current(current_pos(current_round)) = [];
                    remaining_message{current_round} = remaining_message_current;
                end
            end
            
          else
              num_remainings = num_remainings + 1;
              remaining_roots(num_remainings) = current_root; % need to be further investigated, more than one path
          end
        end
        
        
        if (num_remainings==num_remainings_last)
          break;
        end

        num_remainings_last = num_remainings;
        remaining_roots_elimination = current_remaining_roots_elimination;

        % peel-off procedure
        for current_remain_root = remaining_roots
            current_possible_message = survival_paths{current_remain_root, num_round};
            for current_round = 2 : num_round
                current_remaining_messages = remaining_message{current_round};
                current_column = current_possible_message(:, current_round);
                if (isempty(current_column))
                    break;
                end
                num_candidate_current_round = length(current_column);
                num_elimination = 0;
                eliminated_messages = [];
                for l = 1 : num_candidate_current_round
                    if (sum(current_remaining_messages == current_column(l))==0)
                        num_elimination = num_elimination + 1;
                        eliminated_messages(num_elimination) = l;
                    end
                end
                current_possible_message(eliminated_messages, :) = [];
            end
            
            [num_survival_paths(current_remain_root, num_round), ~] = size(current_possible_message);
            survival_paths{current_remain_root, num_round} = current_possible_message;

        end       
   end
   
   flag_resolve_repeated_root = 1;

   % repeated-root case, i.e, N trees with the same root have N survival paths.
   num_duplicate_roots = length(duplicate_decoded_roots);
   for num_group = 1 : num_duplicate_roots
       current_roots = duplicate_decoded_roots{num_group};
       size_current_roots = length(current_roots);
       current_survival_paths = survival_paths{current_roots(1), num_round};
       if (isempty(current_survival_paths))
           continue;
       end
       second_blocks = current_survival_paths(:, 2);
       size_second_blocks = length(second_blocks);
       if (size_second_blocks < 2)
           continue;
       end
       
       [~, unique_second_blocks_idx] = unique(second_blocks, 'stable');
       repeated_sets = [];
       repeated_idx = setdiff(1:size_second_blocks, unique_second_blocks_idx);
       for r_idx = repeated_idx
           delta_set = find(second_blocks==second_blocks(r_idx));
           repeated_sets = [repeated_sets, delta_set'];
       end
       if (~isempty(repeated_sets))
           repeated_sets = unique(repeated_sets, 'stable');
       end
       appear_once_sets = setdiff(1:size_second_blocks, repeated_sets);
       
       if (~isempty(appear_once_sets))
           size_once_appear = length(appear_once_sets);
           flag_resolve_repeated_root = 0;
           idx_r = 0;
           if (size_once_appear > size_current_roots)
               size_once_appear = size_current_roots; % just pick up the first size_current_roots paths.
           end
           
           for root_current = current_roots(1:size_once_appear)
               idx_r = idx_r + 1;
               num_survival_paths(root_current, num_round) = 1;
               survival_paths{root_current, num_round} = current_survival_paths(appear_once_sets(idx_r), :);
           end           
           for root_current = current_roots(size_once_appear+1:size_current_roots)
               num_survival_paths(root_current, num_round) = size_second_blocks-size_once_appear;
               survival_paths{root_current, num_round} = current_survival_paths(repeated_sets, :);
           end           
       end       
   end
   
   if (flag_resolve_repeated_root) % didn't resolve any repeated-root case;
       break;
   end
   
end 

end
function [current_blocks, check] = outer_CRC_checks_BMST_efficient(previous_blocks, previous_superimposed_sum, current_path, current_block, num_protect_section, m, memory, r, CRC_poly)
% the (efficient) CRC check procedure starting with substracting BMST;
% store previous original_block_bits so that the same info doesn't need to compute twice.
% start from the second layer.

% current_blocks: the recovered info bits up to current blocks;
% check: indicator regarding whether the current bits satisfy the CRC condtion;

% previous_blocks: all previous recovered info bits;
% previous_superimposed_sum: the previoused computed bit-wise sum;
% current_path: the current path including all nodes;
% current_block: the current stage or round;
% num_protect_section: the number of protected sections;
% m: the number of bits for each block, i.e., log2(M);
% memory: the memory parameter used in BMST encoding;
% r: a vector indicating the number of redundant bits at each round;
% CRC_poly: the chosen CRC generator polynomial;

    check = 0;
    original_block_bits = zeros(current_block, m);
    original_block_bits(1:current_block-1, :) = previous_blocks;


    current_block_bits = de2bi(current_path(current_block)-1, m);
    original_block_bits(current_block, :) = mod(current_block_bits-previous_superimposed_sum, 2);

    current_path_info_bits = [];
    start_cc = 0;

    for cc = current_block - num_protect_section+1 : current_block-1
      end_cc = start_cc + m - r(cc);
      current_path_info_bits(start_cc+1: end_cc) = original_block_bits(cc, 1:m-r(cc));
      start_cc = end_cc;
    end


    current_blocks = original_block_bits;

    check_bits = [current_path_info_bits, original_block_bits(current_block, :)];
    check = CRC_check(check_bits, CRC_poly);

end



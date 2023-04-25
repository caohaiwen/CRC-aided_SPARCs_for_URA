function [list, num_checks] = outerDecoder(symbol_list, C, parity_lengths, B, LT)
    % Produce a list of messages from a list of symbols and corresponding
    % parity check matrices
    %    
    candidate_paths = symbol_list{1};
    L               = length(parity_lengths);
    num_checks = 0;
    for l = 2:L
        symbols = symbol_list{l};
        new_candidates  = [];
        for i_path = 1:size(candidate_paths,2)
            path = candidate_paths(1:l-1,i_path);
            for i_symbol = 1:length(symbols)
                symbol          = symbols(i_symbol);
                parity_symbol   = rem(symbol,2^parity_lengths(l));
                data_symbol     = bitshift(symbol, -parity_lengths(l));
                p               = parity_check_last([path;data_symbol], parity_symbol,C,log2(B),parity_lengths, LT);
                num_checks = num_checks + 1;
                if p ==0
                    new_candidates = [new_candidates [path;data_symbol]];
                    if size(new_candidates,2) > 4000
                        list = [];
                        disp('Too many paths..');
                        return;
                    end
                end
            end  
        end
        candidate_paths = new_candidates;
        n_paths = size(candidate_paths,2);
        
    end
    list = candidate_paths;
end


% check only the last parity section 
function out = parity_check_last(data_symbols, parity_symbol, C, ~, parity_lengths, LT)
    L               = length(data_symbols);
    parity_lengths  = parity_lengths(1:L);
    out             = 0;
    for i = 1:parity_lengths(L)
        parity_bit = LUT2(bitand(C{L}(i,:)',data_symbols), LT);
        if rem(parity_symbol,2)~= parity_bit
            out = 1;
            break;
        end
        parity_symbol = bitshift(parity_symbol,-1);
    end

end

% Lookuptable for sums of bits
% Input:  list of integer numbers
% Output: 0 if the sum of the bit represented is even
%         1 if its odd
function out = LUT2(symbol_array, LT)
    binary = LT(symbol_array+1);
    out = mod(sum(sum(binary)),2);
end
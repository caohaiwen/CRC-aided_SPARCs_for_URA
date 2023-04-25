function [ sparc_messages, C ] = outerEncoder( message_matrix, B, L, data_lengths, C )
%OUTERENCODER Summary of this function goes here
%   Detailed explanation goes here
    parity_lengths = log2(B) - data_lengths;
    if nargin < 5
        C = createParityMatrices(log2(B),parity_lengths);
    end
    sparc_messages  = encode(message_matrix, C, L, parity_lengths);
    sparc_messages = sparc_messages + 1;
end

function parity_m = createParityMatrices(b,parity_lengths)
    L = length(parity_lengths);
    data_lengths = b - parity_lengths;
    for i = 1:L
        data_bits   = i*b - sum(parity_lengths(1:i)); 
        parity_bits = parity_lengths(i);
        for k = 1:i
            parity_m{i}(:,k) = randi([0,2^data_lengths(k) - 1],parity_bits, 1);
        end
    end
end



function coded_matrix = encode(message_matrix, C, L, parity_lengths)
    K_a             = size(message_matrix,2);
    coded_matrix    = zeros(L, K_a);
    
    for i = 1:L
        data_symbols    = message_matrix(1:i,:);
        parity_symbols  = zeros(1,K_a);
        for k = 1:parity_lengths(i)
            for j = 1:K_a
                parity_bit   = LUT(bitand(C{i}(k,:)',data_symbols(:,j)));
                parity_symbols(j)  = parity_symbols(j) + parity_bit.*2^(parity_lengths(i));
                parity_symbols(j)  = bitshift(parity_symbols(j),-1);
            end

            
        end
        
        coded_matrix(i,:) = bitshift(data_symbols(i,:), parity_lengths(i)) + parity_symbols;
    end
end

% Lookuptable for sums of bits
% Input:  list of integer numbers
% Output: 0 if the sum of the bit represented is even
%         1 if its odd
function out = LUT(symbol_array)
    out = mod(sum(sum(dec2bin(symbol_array) - '0')),2);
end


% covert a matrix where the columns are binary representation of integers
% into integers
function out = bit2dec(bit_matrix)
    b   = size(bit_matrix,1);
    K   = size(bit_matrix,2);
    out = zeros(1,K);
    for i = 1:b
        out = out + bit_matrix(i,:).*2^(i-1);
    end
end
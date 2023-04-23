function c = CRC_encoding(m, g)

% m : information vector
% g : generator polynomial
% c : the CRC part
% since the order in gfdeconv is ascending, we need do an "fliplr"
% operation.

c = zeros(1, length(g)-1);
m_r = [m zeros(1, length(g)-1)];
[~, r] = gfdeconv(fliplr(m_r), fliplr(g));
c(1:length(r)) = r;
c = fliplr(c);

end

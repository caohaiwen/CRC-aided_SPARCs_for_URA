function check = CRC_check(r_c, g)

check = 0;
[~, r] = gfdeconv(fliplr(r_c), fliplr(g));
if (sum(r)==0)
    check = 1;
end

end

function all_rm = mk_all_Pb(encoder,bits)
i = 1;
RM = [];
    while(i <= size(bits,2))
        this_bits = bits(:,i);
        [P, b] = encoder.makePb(this_bits);
        rm = encoder.gen_chirp(P,b.');
        RM = [RM rm];
        i = i + 1;
    end
    
    disp(RM);
    all_rm = RM;
    rev_rm = ctranspose(all_rm);
    mat = all_rm * rev_rm;
    disp(mat);
end
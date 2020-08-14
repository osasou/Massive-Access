
function [propfound, ave_dist] = run(re, m, p, EbN0, K, trials)
    addpath('utils');
%master testing function, runs multiple trials, for definitions of
%parameters see the Encoder and Decoder classes

    sumpropfound=0;
    sum_dist=0;

    for trial = 1:trials
        [propfound_trial, h_dist] = chirrup_test(re,m,p,K,EbN0,[]);
        sumpropfound=sumpropfound+propfound_trial;
        sum_dist=sum_dist+h_dist;
    end

    propfound = sumpropfound/trials;
    ave_dist = sum_dist/trials;
end



function [propfound_trial, h_dist_trial] = chirrup_test(re,m,p,K,EbN0,input_bits)
    %runs a single test
    %;‚ğ‚Â‚¯‚é‚Æˆ—‚ğ‰B‚¹‚é
    h_all_in = [];
    h_all_out = [];
    encoder=Encoder(re,m,p,K,EbN0,input_bits);
    [encoder, input_bits] = encoder.generate_random_bits(); %we update the value of input_bits to save confusion, despite already being encoder.input_bits
    [Y, h_all_in] = encoder.chirrup_encode;
    decoder=Decoder(Y,re,m,p,K);
    [output_bits, h_all_out] = decoder.chirrup_decode();
    propfound_trial = compare_bits(input_bits,output_bits);
    h_dist_trial = square_dist(h_all_in, h_all_out);
end



encoder = Encoder(0,5,0,1,200,NaN)
A = [0:2^encoder.B-1]
b = de2bi(A);
bits = flip(b);
bits = flip(b, 2)
bits = bits.';
mk_all_Pb(encoder,bits)

とりあえずMassive_accessにあるアルゴリズムで作成してみたが、うまくはいかないみたい
ユーザが一人の時は10回に1回くらいでうまくいく
雑音なしでやると途中でyの値が0になってしまって上手くいかない。

とりあえずここで0_1は終わり


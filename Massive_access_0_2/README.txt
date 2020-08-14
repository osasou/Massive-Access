

コマンドウィンドウ上で
encoder = Encoder(0,3,0,1,200,NaN)
A = [0:2^encoder.B-1]
b = de2bi(A);
bits = flip(b);
bits = flip(b, 2)
bits = bits.';
[P, b] = encoder.makePb(bits(:,1))
rm = encoder.gen_chirp(P,b.')
とやれば、RM符号の符号語がわかる

とりあえずMassive_accessにあるアルゴリズムで作成してみたが、うまくはいかないみたい
ユーザが一人の時は10回に1回くらいでうまくいく
雑音なしでやると途中でyの値が0になってしまって上手くいかない。

とりあえずここで0_1は終わり

次 0_2 

Pee1 =

     1     0     0
     0     1     0
     0     0     1

bee1 =

   1   0   0

の時には上手くいかない

Pee1 =

     1     0     1
     0     1     1
     1     1     0

bee1 =

   1   1   0

はPはちゃんと復号できたけど b は 0 1 0 となる

bit = [1 0 1 1 1 0 1 1 0]
[P, b] = encoder.makePb(bit)
rm = encoder.gen_chirp(P,b)

decoder = Decoder({rm},0,3,0,1)
decoder.chirrup_decode




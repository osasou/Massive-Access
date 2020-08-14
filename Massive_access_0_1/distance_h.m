function x = distance_h(in, out)
   %%inがslotを分けた後だから、3ユーザの場合6つのhの値。outがslotを含め、複合した後だから、3つのhの値。
   %%→じゃあ、Decodeの時に6つの時はない？　　わからない。。。
   
   x = in - out;
%    x = 0;
end
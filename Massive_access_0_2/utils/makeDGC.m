function K = makeDGC(re,m)

[i j] = upper_indices(m,re);
%disp([i,j])
%makeDGC(1,3)の場合、
%i=[1 1 2]
%j=[2 3 3]と入る

for p = 1:length(i)
    
    K(:,:,p) = zeros(m,m);
    % i行j列とj行i列に1を代入する
    % K_1には、1行2列と2行1列
    % K_2には、1行3列と3行1列
    % K_3には、2行3列と3行2列
    K(i(p),j(p),p) = 1;
    K(j(p),i(p),p) = 1;
    
end

%disp(K)
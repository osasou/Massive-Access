function dist = square_dist(h_all,h_hat)

sum = 0;

% for i=1:length(h_all)
%     for j=1:length(h_hat)
%         % sum = sum + (h_all(i) - h_hat(j))^2;
%         % ↑だと{(a+bi)*(c+bi)}^2になる
%     end
% end

% sum = sum/length(h_all);

% err = sqrt(sum);
sum = h_all - h_hat;
dist = abs(sum);
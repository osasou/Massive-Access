function [x,y] = random_circle(range_max,num)
%RANDOM_CIRCLE この関数の概要をここに記述
%   詳細説明をここに記述


% rng('shuffle')
x = [];
y = [];
a = -range_max;
b = range_max;
cnt = 1;
while (cnt <= num)
    num_x = (b-a).*rand() + a;
    num_y = (b-a).*rand() + a;

    if num_x^2+num_y^2 <= range_max^2
        
        x = [x,num_x];
        y = [y,num_y];
        %{
        x = [x,sqrt(2)/2];
        y = [y,sqrt(2)/2];
        %}
        cnt = cnt+1;
    end
end


%scatter(x,y);

%{
rng(0,'twister')
rvals = 2*rand(1000,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(1000,1);
radii = 100*(rand(1000,1).^(1/3));
[x,y,z] = sph2cart(azimuth,elevation,radii);
figure
plot3(x,y,z,'.')
axis equal
%}
end


function [ Y_min, beta_c, X_min, contLoc  ] = min_point_ellipse( pos, ori, ellipse, res )
%MIN_POINT_ELLIPSE Summary of this function goes here
%   Detailed explanation goes here

x = pos(1);
y = pos(2);

theta = ori;

a_e = ellipse(1);
b_e = ellipse(2);

beta = linspace(0,2*pi,res)';

r = a_e*b_e./(sqrt(b_e^2*cos(beta).*cos(beta)+a_e^2*sin(beta).*sin(beta)));

v = r .* cos(beta);
w = r .* sin(beta);

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

ic = zeros(2,length(beta));

for i=1:length(beta)
    
    ic(:,i) = [x;y] + R*[v(i);w(i)];
    
end

[Y_min, ind ] = min(ic(2,:));
beta_c = beta(ind);

% do fminunc version, starting from discretely computed beta_c 
options = optimset('display', 'off', 'TolFun', 1e-8, 'LargeScale', 'off');
beta_c = fminunc(@(x) minimum_point_f(x,ori,ellipse), beta_c, options);

% scale u as necessary
r = a_e*b_e./(sqrt(b_e^2*cos(beta_c).*cos(beta_c)+a_e^2*sin(beta_c).*sin(beta_c)));
u = [r*cos(beta_c) r*sin(beta_c)]'; 
Y_min = pos(2) + R(2,:)*u;
X_min = pos(1) + R(1,:)*u;
contLoc = [R(1,:)*u,R(2,:)*u];
end


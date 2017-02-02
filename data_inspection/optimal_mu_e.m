function [ mu_opt,e_opt,P_estimate,P_measured,l_pre ] = optimal_mu_e( model_str, trial, with_text )
%OPTIMAL_MU_E Summary of this function goes here
%   Detailed explanation goes here

addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/models')
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/data')
load('ellipse_uniform')

m = 36.4*1e-3; g = 9.81;

% Define size of square
a_e = 70/1000/2;
b_e = 50/1000/2;
I_inertia = m*(b_e^2+a_e^2)/5;
Mass = [m,0,0;
    0,m,0;
    0,0,I_inertia];

v = bounce_array(trial).states(4:6)';
d = bounce_array(trial).d;
n = bounce_array(trial).n;
h = 1/250;
ha = h*[0;-g;0];

vm = bounce_array(trial).states(10:12)';


J = [d;n];
Minv = J*(Mass\J');
M    = inv(Minv);

rad = (J*v)'*M*(J*v);
MVI = M*(J*v);

Pt = M*((J*vm)-(J*v));
vpost = M*(J*vm);

% mu_init = [0.01,0.1,0.1,0.1,0.2,0.3];
% e_init  = [0.60,0.6,0.7,0.9,0.6,0.6];

mu_init = [0.01,0.1,0.2,0.3];
e_init  = [0.60,0.7,0.6,0.6];

% mu_init = 0.1;
% e_init  = 0.6;

mu_lim = [0,0.4];
ep_lim = [0.4,1];
lb = [mu_lim(1),ep_lim(1)];
ub = [mu_lim(2),ep_lim(2)];

% options = optimoptions('fmincon','FiniteDifferenceType','central',...
%     'StepTolerance',1e-10,'Display','iter-detailed');
options = optimoptions('fmincon','FiniteDifferenceType','central',...
    'StepTolerance',1e-10,'Display','off');
opt_vec = zeros(size(mu_init,2),3);

for opt_count = 1:size(mu_init,2)
    tic
    x0 = [mu_init(opt_count),e_init(opt_count)];
    [x, fval] = fmincon(@(x)cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt), x0, ...
        [], [], [], [], lb, ub, [], options);
    opt_vec(opt_count,:) = [x(1),x(2),fval];
    
    if with_text
        fprintf('Attempt %d:\n',opt_count)
        fprintf('mu_opt: %2.2f  eps_opt: %2.2f   fval: %2.4f   opt_time: %2.4f \n',x(1),x(2),fval,toc)
    end
end

[~,ind]=min(opt_vec(:,3));
mu_opt = opt_vec(ind,1);
e_opt  = opt_vec(ind,2);


[ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu_opt, e_opt );
p_hat = (M*(J*(v1-v)))';

P_estimate = p_hat;
P_measured = Pt;
l_pre = (J*v);
end

function obj = cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt)

mu = x(1);
e  = x(2);
[ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
p_hat = (M*(J*(v1-v)))';
obj = norm(Pt-p_hat');

end
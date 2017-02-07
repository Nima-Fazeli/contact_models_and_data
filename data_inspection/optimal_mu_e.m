function [ mu_opt,e_opt,P_estimate,P_measured,l_pre ] = optimal_mu_e( model_str, trial, opts )
%OPTIMAL_MU_E Summary of this function goes here
%   Detailed explanation goes here
with_text   = opts{1};
solver_name = opts{2};

addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/models')
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/data')
load('ellipse_uniform')

if sum(bounce_array(trial).flags)<1

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

mu_init = [0.01,0.1,0.1,0.1,0.2,0.3];
e_init  = [0.60,0.6,0.7,0.9,0.6,0.6];

% mu_init = [0.01,0.1,0.2,0.3];
% e_init  = [0.60,0.7,0.6,0.6];

% mu_init = 0.1;
% e_init  = 0.6;

mu_lim = [0,0.4];
ep_lim = [0.2,1];
lb = [mu_lim(1),ep_lim(1)];
ub = [mu_lim(2),ep_lim(2)];

switch solver_name
    case 'GA'
        tic
        nvars = 2;
        [x,fval] = ga(@(x)cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt),nvars,[],...
            [],[],[],lb,ub,[]);
        mu_opt = x(1);
        e_opt  = x(2);
        if with_text
            fprintf('mu_opt: %2.2f  eps_opt: %2.2f   fval: %2.4f   opt_time: %2.4f \n',x(1),x(2),fval,toc)
        end
        
    case 'fmincon'

        options = optimoptions('fmincon','FiniteDifferenceType','central',...
            'StepTolerance',1e-10,'Display','off');
        tic
        x0 = opts{3};
        [x, fval] = fmincon(@(x)cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt), x0, ...
            [], [], [], [], lb, ub, [], options);
        
        if with_text
            fprintf('mu_opt: %2.2f  eps_opt: %2.2f   fval: %2.4f   opt_time: %2.4f \n',x(1),x(2),fval,toc)
        end
        
        mu_opt = x(1);
        e_opt  = x(2);
        
    case 'global_search'
        opt_options = optimoptions(@fmincon,'Algorithm','interior-point');
        tic
        problem = createOptimProblem('fmincon','objective',...
            @(x)cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt),'x0',[0.1,0.6],...
            'lb',lb,'ub',ub,'options',opt_options);
        gs = GlobalSearch;
        [x,fval,~,~,solutions]=run(gs,problem);
        mu_opt = solutions(1).X(1);
        e_opt  = solutions(1).X(2);
        if with_text
            fprintf('mu_opt: %2.2f  eps_opt: %2.2f   fval: %2.4f   opt_time: %2.4f \n',x(1),x(2),fval,toc)
        end
        
end

[ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu_opt, e_opt );
p_hat = (M*(J*(v1-v)))';

P_estimate = p_hat;
P_measured = Pt;
l_pre = (J*v);

else
    mu_opt=-1;
    e_opt =-1;
    P_estimate=zeros(2,1);
    P_measured=zeros(2,1);
    l_pre =zeros(2,1);
end

end

function obj = cost_function(x,model_str,Mass,M,n,d,v,ha,J,Pt)

mu = x(1);
e  = x(2);
[ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
p_hat = (M*(J*(v1-v)))';
obj = norm(Pt-p_hat');

end
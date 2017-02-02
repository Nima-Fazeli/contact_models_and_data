function [  ] = convexity_check( model_str, trial )

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


e_range  = linspace(0.40,0.9,20);
mu_range = linspace(0.01,0.4,20);
[E,Mu] = meshgrid(e_range,mu_range);
norm_error = zeros(size(E,1),size(E,2));

p_list  = zeros(2,size(E,1)*size(Mu,1));
counter = 0;
for ii=1:size(E,1)
    for jj=1:size(Mu,1)
        counter = counter+1;
        mu = Mu(ii,jj);
        e  = E(ii,jj);
        [ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
        p_hat = (M*(J*(v1-v)))';
        p_list(:,counter) = p_hat;
        norm_error(ii,jj) = sqrt(sum((Pt-p_hat').^2,1));
    end
end

figure
surf(E,Mu,norm_error)
xlabel 'Epsilon' 
ylabel 'Mu' 
zlabel '2 Norm Error'

figure
hold on
plot(p_list(1,:),p_list(2,:),'ob')
plot(Pt(1),Pt(2),'or')

% vec_test = sqrt(sum((p_list-Pt).^2,1));
% figure
% plot(vec_test)

end


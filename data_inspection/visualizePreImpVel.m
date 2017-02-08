function [  ] = visualizePreImpVel(  )

close all
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/models')
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/data')
load('ellipse_uniform')

%%% only using this to plot results for data that is clean with no flags
load('epsilon_vy_APNewton')
ind_to_plot = find(data_vec(:,1)>-0.5);

%%% data preperation
V=cellfun(@(x) x(4:6),{bounce_array(1:2000).states},'UniformOutput',0);
N=cellfun(@(x) x(1:3),{bounce_array(1:2000).n},'UniformOutput',0);
D=cellfun(@(x) x(1:3),{bounce_array(1:2000).d},'UniformOutput',0);

V = cell2mat(V);
V = reshape(V,3,[]);
N = cell2mat(N);
N = reshape(N,3,[]);
D = cell2mat(D);
D = reshape(D,3,[]);

Vn = dot(V,N);
Vt = dot(V,D);

%%% plot script
figure
plot(Vt(ind_to_plot),Vn(ind_to_plot),'o')
xlabel('Tangential pre-impact velocity of contact point')
ylabel('Normal pre-impact velocity of contact point')
grid on
end


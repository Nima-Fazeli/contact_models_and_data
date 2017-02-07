function [  ] = visualizeMuEpsilon(  )

close all
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/models')
addpath('/Users/nimafazeli/Documents/MATLAB/2017 01 - Learning Contact New Data/data')
load('epsilon_vy_APNewton')
% load('epsilon_vy_APPoisson')
% load('epsilon_vy_mirtich')
% load('epsilon_vy_wang')
% load('epsilon_vy_DrumShell')

ind_to_plot = find(data_vec(:,1)>-0.5);
figure
subplot(211)
plot(-data_vec(ind_to_plot,3),data_vec(ind_to_plot,1),'o')
xlabel 'Linear Velocity pre-impact X'
ylabel 'Coeff of Friction'
subplot(212)
plot(-data_vec(ind_to_plot,4),data_vec(ind_to_plot,2),'o')
xlabel 'Linear velocity pre-impact Y'
ylabel 'Coeff of Restitution'

figure
subplot(211)
plot(-data_vec(ind_to_plot,5),data_vec(ind_to_plot,1),'o')
xlabel 'Change in Linear Momentum X'
ylabel 'Coeff of Friction'
subplot(212)
plot(-data_vec(ind_to_plot,6),data_vec(ind_to_plot,2),'o')
xlabel 'Change in Linear Momentum Y'
ylabel 'Coeff of Restitution'

figure
subplot(211)
ksdensity(data_vec(ind_to_plot,1))
ylabel('PDE of Mu')
xlabel('Mu Values')
grid on
subplot(212)
ksdensity(data_vec(ind_to_plot,2))
ylabel('PDE of Epsilon')
xlabel('Epsilon Values')
grid on
end


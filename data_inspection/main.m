
clc
clear
close all

mdl_id   = 3;
trial_id = 3;

mdl_list={'DrumShell','APPoisson','APNewton','mirtich','wang','whittaker'};
fprintf('Looking at results from model %s on trial number %d \n',mdl_list{mdl_id},trial_id)

% %%% visualize the energy ellipse and hypotheses
% visualizeEnergyEllipse(mdl_list{mdl_id},trial_id)
% 
% %%% visualize the convexity of of the error in momentum
% convexity_check(mdl_list{mdl_id},trial_id)
% 
% %%% estimate optimal mu and epsilon based on minimizing contact point
% %%% momentum error
% with_text = 1;
% [mu,e,~,~,~]=optimal_mu_e(mdl_list{mdl_id},trial_id, with_text);

data_num = 50;
data_vec = zeros(data_num,4);
with_text = 0;
for i=1:data_num
    [mu,e,~,~,l_pre]=optimal_mu_e(mdl_list{mdl_id}, i, with_text);
    data_vec(i,:) = [mu,e,l_pre'];
end

figure
plot(data_vec(:,4),data_vec(:,2),'o')
xlabel 'Linear Momentum pre-impact'
ylabel 'Coeff of Res.'
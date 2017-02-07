
clc
clear
close all

mdl_id   = 4; % picks the model to be tested
trial_id = 2; % picks the trial to look at.

mdl_list={'DrumShell','APPoisson','APNewton','mirtich','wang','whittaker'};
fprintf('Looking at results from model %s on trial number %d \n',mdl_list{mdl_id},trial_id)

% %%% visualize the energy ellipse and hypotheses
visualizeEnergyEllipse(mdl_list{mdl_id},trial_id)
 
% %%% visualize the convexity of of the error in momentum
with_vis = 1;
ic_mue     = convexity_check(mdl_list{mdl_id},trial_id,with_vis);
ic_mue_com = convexity_check_COM(mdl_list{mdl_id},trial_id,with_vis);
% %%% estimate optimal mu and epsilon based on minimizing contact point
% %%% momentum error
with_text = 1;
options_mue = {with_text,'fmincon',ic_mue};
[mu,e,~,~,~]=optimal_mu_e(mdl_list{mdl_id},trial_id, options_mue);

%%
data_num = 500;
data_vec = zeros(data_num,6);
with_text = 0;
options_mue = {with_text,'fmincon',ic_mue};
with_vis = 0;
for i=1:data_num
    ic_mue = convexity_check(mdl_list{mdl_id},i,with_vis);
    options_mue = {with_text,'fmincon',ic_mue};
    [mu,e,~,Pt,l_pre]=optimal_mu_e(mdl_list{mdl_id}, i, options_mue);
    data_vec(i,:) = [mu,e,l_pre',Pt'];
    
    if rem(i,10)==0
        fprintf('Percent done ... %2.0f\n',i/data_num*100);
    end
end

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
xlabel 'Linear Momentum X'
ylabel 'Coeff of Friction'
subplot(212)
plot(-data_vec(ind_to_plot,6),data_vec(ind_to_plot,2),'o')
xlabel 'Linear Momentum Y'
ylabel 'Coeff of Restitution'
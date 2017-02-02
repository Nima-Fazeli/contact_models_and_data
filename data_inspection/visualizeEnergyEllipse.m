function [  ] = visualizeEnergyEllipse( model_str, trial )
%VISUALIZEENERGYELLIPSE Summary of this function goes here
%   Detailed explanation goes here

close all
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

v_list = zeros(50,3);
P_list_drum = zeros(50,2);
P_list_whit = zeros(50,2);
P_list_wang = zeros(50,2);
P_list_appo = zeros(50,2);

mu_range = [0.001,0.3];
ep_range = [0.500,0.8];

for ii=1:500
    mu =  (mu_range(2)-mu_range(1))*rand(1)+ mu_range(1);
    e  =  (ep_range(2)-ep_range(1))*rand(1)+ ep_range(1);
    [ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
    v_list(ii,:) = v1';
    P_list_drum(ii,:) = (M*(J*(v1-v)))';
end

%%% veloctiy
% [xe,ye,~] = drawE2(Minv, zeros(2,1), rad);
% mu = 0.35;
% figure
% hold on
% plot(xe,ye,'b')
% [~,indt] = min(xe);
% [~,indp] = max(xe);
% [~,indtt] = min(ye);
% [~,indpp] = max(ye);
% xaxp=linspace(0,max(xe),100);
% xaxn=linspace(min(xe),0,100);
% plot(xaxp,1/mu*xaxp,'r')
% plot(xaxn,-1/mu*xaxn,'r')
% plot(vpost(1),vpost(2),'ok')
% plot([xe(indt),xe(indp)],[ye(indt),ye(indp)],'--k')
% plot([xe(indtt),xe(indpp)],[ye(indtt),ye(indpp)],'--k')
% for ii=1:50
% plot(v_list(ii,1),v_list(ii,2),'or')
% % plot(P_list_whit(ii,1),P_list_whit(ii,2),'*b')
% % plot(P_list_wang(ii,1),P_list_wang(ii,2),'^r')
% % plot(P_list_appo(ii,1),P_list_appo(ii,2),'^g')
% end
% xlabel 'CP V x'
% ylabel 'CP V y'
% box
% grid
% axis equal

%%% impulse
[xe,ye,~] = drawE2(Minv, MVI, rad);
mu = 0.35;
figure
hold on
plot(xe,ye,'b')
[~,indt] = min(xe);
[~,indp] = max(xe);
[~,indtt] = min(ye);
[~,indpp] = max(ye);
xaxp=linspace(0,max(xe),100);
xaxn=linspace(min(xe),0,100);
plot(xaxp,1/mu*xaxp,'r')
plot(xaxn,-1/mu*xaxn,'r')
plot(Pt(1),Pt(2),'ok')
plot([xe(indt),xe(indp)],[ye(indt),ye(indp)],'--k')
plot([xe(indtt),xe(indpp)],[ye(indtt),ye(indpp)],'--k')
for ii=1:50
plot(P_list_drum(ii,1),P_list_drum(ii,2),'or')
% plot(P_list_whit(ii,1),P_list_whit(ii,2),'*b')
% plot(P_list_wang(ii,1),P_list_wang(ii,2),'^r')
% plot(P_list_appo(ii,1),P_list_appo(ii,2),'^g')
end
xlabel 'CP Momentum x'
ylabel 'CP Momentum y'
box
grid
axis equal


%%% keep epsilon constant

e = 0.7;
range = linspace(-6,0,100);
mu_list = exp(range);
P_list_drum =zeros(length(mu_list),2);
v_list_drum =zeros(length(mu_list),3);

for i=1:length(mu_list)
    mu = mu_list(i);
    [ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
    P_list_drum(i,:) = (M*(J*(v1-v)))';
    v_list_drum(i,:) = v1';
end

figure
hold on
plot(xe,ye,'b')
[~,indt] = min(xe);
[~,indp] = max(xe);
xaxp=linspace(0,max(xe),100);
xaxn=linspace(min(xe),0,100);
plot(xaxp,1/mu*xaxp,'r')
plot(xaxn,-1/mu*xaxn,'r')
plot(Pt(1),Pt(2),'ok')
plot([xe(indt),xe(indp)],[ye(indt),ye(indp)],'--k')
for ii=1:length(mu_list)
    xp = P_list_drum(ii,1);
    yp = P_list_drum(ii,2);
    plot(xp,yp,'or')
    txt1 = ['\leftarrow mu =',num2str(mu_list(ii)),''];
%     text(xp,yp,txt1)
end
xlabel 'CP Momentum x'
ylabel 'CP Momentum y'
box
grid
axis equal
title 'Mu variation'

figure
subplot(311)
plot(mu_list)
ylabel 'mu values tried'
subplot(312)
plot(mu_list,v_list_drum(:,1));
ylabel 'Tangent velocity versus mu'
subplot(313)
plot(mu_list,P_list_drum(:,1))
ylabel 'Tangent momentum versus mu'

figure
subplot(121)
plot(v_list_drum(:,1),v_list_drum(:,2),'o')
ylabel 'CP Velocity y'
xlabel 'CP Velocity x'
title 'Mu variation'
subplot(122)
plot(v_list_drum(:,1)-v(1),v_list_drum(:,2)-v(2),'o')
ylabel 'Change CP Velocity y'
xlabel 'Change CP Velocity x'
title 'Mu variation'

%%% keep mu constant

mu = 0.2;
e_list = linspace(0.5,0.9,20);
P_list_drum =zeros(length(e_list),2);
v_list_drum =zeros(length(e_list),3);

for i=1:length(e_list)
    e = e_list(i);
    [ v1, ~ ] = model_eval( model_str, Mass, n, d, v, ha, mu, e );
    P_list_drum(i,:) = (M*(J*(v1-v)))';
    v_list_drum(i,:) = v1';
end

figure
hold on
plot(xe,ye,'b')
[~,indt] = min(xe);
[~,indp] = max(xe);
xaxp=linspace(0,max(xe),100);
xaxn=linspace(min(xe),0,100);
plot(xaxp,1/mu*xaxp,'r')
plot(xaxn,-1/mu*xaxn,'r')
plot(Pt(1),Pt(2),'ok')
plot([xe(indt),xe(indp)],[ye(indt),ye(indp)],'--k')
for ii=1:length(e_list)
    xp = P_list_drum(ii,1);
    yp = P_list_drum(ii,2);
    plot(xp,yp,'or')
    txt1 = ['\leftarrow e =',num2str(e_list(ii)),''];
%     text(xp,yp,txt1)
end
xlabel 'CP Momentum x'
ylabel 'CP Momentum y'
box
grid
axis equal
title 'Epsilon variation'

figure
subplot(311)
plot(e_list)
ylabel 'epsilon values tried'
subplot(312)
plot(e_list,v_list_drum(:,1));
ylabel 'Tangent velocity versus epsilon'
subplot(313)
plot(e_list,P_list_drum(:,1))
ylabel 'Tangent momentum versus epsilon'

figure
subplot(121)
plot(v_list_drum(:,1),v_list_drum(:,2),'o')
ylabel 'CP Velocity y'
xlabel 'CP Velocity x'
title 'Epsilon variation'
subplot(122)
plot(v_list_drum(:,1)-v(1),v_list_drum(:,2)-v(2),'o')
ylabel 'Change CP Velocity y'
xlabel 'Change CP Velocity x'
title 'Epsilon variation'

end


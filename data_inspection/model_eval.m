function [ v, z ] = model_eval( model_str, Mass, n, d, v, ha, mu, e )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp(model_str,'DrumShell')
    [v, z] = DrumShell(Mass, n', d', v, ha, mu, e);
end

if strcmp(model_str,'whittaker')
    [v, z] = whittaker(Mass, n', d', v, ha, mu, e);
end

if strcmp(model_str,'wang')
    [v, z] = wang(Mass, n', d', v, ha, mu, e);
end

if strcmp(model_str,'APPoisson')
    [v, z] = APPoisson(Mass, n', d', v, ha, mu, e);
end

if strcmp(model_str,'APNewton')
    [v, z] = APNewton(Mass, n', d', v, ha, mu, e);
end

if strcmp(model_str,'mirtich')
    [v, z] = mirtich(Mass, n', d', v, ha, mu, e);
end

% switch model_str
%     case 'DrumShell'
%         [v, z] = DrumShell(Mass, n', d', v, ha, mu, e);
%     case 'whittaker'
%         [v1, z] = whittaker(Mass, n', d', v, ha, mu, e);
%     case 'wang'
%         [v1, z] = wang(Mass, n', d', v, ha, mu, e);
%     case 'APPoisson'
%         [v1, z] = APPoisson(Mass, n', d', v, ha, mu, e);
%     case 'APNewton'
%         [v1, z] = APNewton(Mass, n', d', v, ha, mu, e);
%     case 'mirtich'
%         [v1, z] = mirtich(Mass, n', d', v, ha, mu, e);
% end


end

%     
%     [v1, ~] = whittaker(Mass, n', d', v, ha, mu, e);
%     P_list_whit(ii,:) = (M*((J*v1)-(J*v)))';
%     
%     [v1, ~] = wang(Mass, n', d', v, ha, mu, e);
%     P_list_wang(ii,:) = (M*((J*v1)-(J*v)))';
%     
%     [v1, ~] = APPoisson(Mass, n', d', v, ha, mu, e);
%     P_list_appo(ii,:) = (M*((J*v1)-(J*v)))';
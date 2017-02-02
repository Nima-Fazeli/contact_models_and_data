function [ x,y,minmax ] = drawE2( m_c, offset, p )
%DRAWE2 Summary of this function goes here
%   Detailed explanation goes here

[v,d] = eig(m_c/p);
a = sqrt(1/d(1,1));
b = sqrt(1/d(2,2));
phi = atan2(v(2,1),v(1,1));

t=linspace(0,2*pi,100);
x=-offset(1)+a*cos(t)*cos(phi)-b*sin(t)*sin(phi);
y=-offset(2)+a*cos(t)*sin(phi)+b*sin(t)*cos(phi);

[~,minpt] = min(x);
[~,maxpt] = max(x);
minmax = [minpt;maxpt];
end


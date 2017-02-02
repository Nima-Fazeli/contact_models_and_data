% function for evaluating the minimum point of an ellipse, for use with
% unconstrained minimimization (fminunc)
% x: the 2d position of the c.o.m. of the ellipse
% theta: the planar orientation of the ellipse
% ellipse: the radii of the ellipse
function y = minimum_point_f(x, theta, ellipse)

  a_e = ellipse(1);
  b_e = ellipse(2);

  % determine u
  S = [cos(x) -sin(x); sin(x) cos(x)];
  u = S*[1 0]';

  % scale u as necessary
  c = sqrt((u(1)^2)/(a_e^2) + (u(2)^2)/(b_e^2));
  u = u/c;

  % form R
  R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

  % get the point
  p = R*u;
  y = p(2);


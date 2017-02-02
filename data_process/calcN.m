% computes the 'n' vector used in Anitescu-Potra for the ellipse
% NOTE: this can be made more efficient
function n = calcN(p, q)
  r = [p(1); p(2); 0] - [q(1); q(2); 0];
  comp = cross(r, [0 1 0]');
  n = [0 1 comp(3)]';


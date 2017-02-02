% computes the 'd' vector used in Anitescu-Potra for the ellipse
% NOTE: this can be made more efficient
function d = calcD(p, q)

  r = [p(1); p(2); 0] - [q(1); q(2); 0];
  comp = cross(r, [1 0 0]');
  d = [1 0 comp(3)]';


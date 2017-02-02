
% u is the contact-space velocity ([tan1 tan2 normal])
% mu is the coefficient of friction
% epsilon is the coefficient of restitution
function [v_plus, z] = whittaker(M, n, s, v, ha, mu, epsilon)

  % get the normal velocity
  nv = n'*(v + ha);

  % if nv is non-negative, quit
  if (nv >= 0)
    z = [0 0 s'*(v + ha)];
    v_plus = v;
    return;
  end

  % get the normal and tangent velocity
  nsv = [n*(1+epsilon) s]'*(v + ha);

  % compute the contact space inertia matrix
  K = [n s]'*(M \ [n s]);

  % compute the full friction impulse
  z = K \ -nsv;

  % friction cone check
  if (abs(z(2)) > mu*abs(z(1)))
    % impulse is outside the friction cone; limit it
    z(2) = mu*z(1);

    % compute the new velocity
    v_plus = M \ [n s]*z + (v + ha);

    % if the tangent velocity and the impulse have the same sign, reverse
    % the tangent impulse
    if (z(2)*(s'*v_plus) > 0)
      z(2) = -z(2);
    end
  end

  % compute the new velocity
  v_plus = M \ [n s]*z + (v + ha);
  z(4) = abs(s'* v_plus);


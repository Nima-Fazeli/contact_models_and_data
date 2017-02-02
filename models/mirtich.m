
% u is the contact-space velocity ([tan1 tan2 normal])
% mu is the coefficient of friction
% epsilon is the coefficient of restitution
function [v_plus, z] = mirtich(M, n, s, v, ha, mu, epsilon)

  % determine K0, u0
  K = [s'*(M \ s), s' * (M \ n);
        n'*(M \ s), n' * (M \ n)];
  u0 = [s'*(v + ha); n'*(v + ha)];
  u = u0;
  
  % setup the step size for integration
  dt = 1e-4; 

  % if the relative velocity is non-negative at the collision point, return w/o applying impulse -- bodies are not impacting
  % NOTE: u(2) is relative normal velocity
  if (u0(2) >= 0)
    z = zeros(3,1);
    v_plus = [M \ s, M \ n]*z + (v + ha);
    z(3) = 0;
    if (z(2) < 0)
      z(3) = -z(2);
      z(2) = 0;
    end
    z(4) = abs(s'*v_plus);
    return;
  end

  % invert K 
  Kinv = pinv(K);
  
  % copy u
  ux = u(1);
  uy = 0;
  uz = u(2);

  % setup work in compression and decompression phases
  Wd = 0;
  Wc = 0;

  % indicate this is the first time running the loop
  firsttime = 1;

  % integration loop
  while (uz < 0 || Wd < -epsilon^2*Wc)

    % check whether ux, uy = 0 (sticking has occurred)
    if (norm([ux uy]) < sqrt(eps))

      % check whether sticking is stable
      if (Kinv(1,2)  <= mu*Kinv(2,2))

        % sticking is stable
        fprintf(1, 'Stable sticking detected\n');

        if (uz < 0)

          % in a compression phase
          Wc = Wc - Kinv(2,2)*uz^2/2;
          uz = 0;

        end

        % in a decompression phase
        uz = sqrt(2/Kinv(2,2)*(-epsilon^2*Wc - Wd) + uz^2);
        ux = 0;

        % compute j and return
        z = K \ ([ux uz]' - u0);
        v_plus = [M \ s, M \ n]*z + (v + ha);
        z(3) = 0;
        if (z(2) < 0)
          z(3) = -z(2);
          z(2) = 0;
        end
        z(4) = abs(s'*v_plus);
        return;

      else

        % sticking is unstable: compute mu
        fprintf(1, 'Unstable sticking detected\n');
        theta = calctheta(K, mu);

        % compute k
        k = K * [-mu; 1];

        % store current uz
        uold = uz;
        if (uz < 0)
          % in a compression phase
          Wc = Wc -uz^2/(2*k(2));
          uz = 0;
        end

        % in a decompression phase
        uz = sqrt(2*k(2)*(-epsilon^2*Wc - Wd) + uz^2);

        % update ux and uy
        ux = ux + k(1)/k(2)*(uz - uold);

        % compute j
        z = K \ ([ux uz]' - u0);
        v_plus = [M \ s, M \ n]*z + (v + ha);
        z(3) = 0;
        if (z(2) < 0)
          z(3) = -z(2);
          z(2) = 0;
        end
        z(4) = abs(s'*v_plus);
        return;

      end
    else % end sticking friction check

      % integrate work
      if (uz < 0)
        Wc = Wc + uz*dt;
      else
        Wd = Wd + uz*dt;
      end

      % update u
      u = u + dt*K*[-mu*ux; 1]; 
      ux = u(1);
      uz = u(2);
    end
  end
 
  % finally, compute the desired impulse -- see [Mirtich, 1996], p. 67
  z = K \ (u - u0);
  v_plus = [M \ s, M \ n]*z + (v + ha);
  z(3) = 0;
  if (z(2) < 0)
    z(3) = -z(2);
    z(2) = 0; 
 end
 z(4) = abs(s'*v_plus);


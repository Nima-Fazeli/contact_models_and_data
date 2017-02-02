% the Anitescu-Potra model with Poisson restitution
%
% @article{Anitescu:1997,
%	Author = {Miha Anitescu and Florian A. Potra},
%	Journal = {Nonlinear Dynamics},
%	Pages = {231--247},
%	Title = {Formulating Dynamic Multi-Rigid-Body Contact Problems with Friction as Solvable Linear Complementarity Problems},
%	Volume = {14},
%	Year = {1997}}

function [v_plus, z] = APPoisson(M, n, s, v, ha, mu, e)

   E = [1;1];
   D = [-s, s];
   
   % Define LCP:
   M_hat = [n'*(M \ n), n'* (M \ D),  0;
       D'*(M \ n), D'*(M \ D),  E;
       mu,       -E',  0];

   q_hat = [n'*(v + ha); D'*(v + ha); 0];
   
   z = LCP(M_hat,q_hat);

   % now add in normal restitution
   z(1) = z(1) * (1 + e);

   % make sure that solution is still ok
   v_plus = [M \ n, M \ D]*[z(1:3)]+(v + ha);
   if (n'*v_plus < 0)
     % re-solve using new velocity
     q_hat = [n'*v_plus; D'*v_plus; 0];
     z = LCP(M_hat, q_hat);
   end
   
   %         if abs(z(1))<0.05
   %             e = 0;
   % %             vlplus1 = zeros(3,1)+1e-3;
   %         end
   


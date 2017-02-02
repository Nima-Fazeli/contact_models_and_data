% the Drumwright-Shell model (Poisson restitution)
% @incollection{Drumwright:2011a,
%	Author = {Evan Drumwright and Dylan A. Shell},
%	Booktitle = {Algorithmic Foundations of Robotics {IX}},
%	Publisher = {Springer Berlin / Heidelberg},
%	Title = {Modeling Contact Friction and Joint Friction in Dynamic Robotic Simulation Using the Principle of Maximum Dissipation},
%	Year = {2011}}
% 
function [v_plus, z] = DrumShell(M, n, s, v, ha, mu, e)

   E = [1;1];
   D = [-s, s];
   
   % Define QP as an LCP: first define quadratic objective matrix
   G = [n'*(M \ n), n'* (M \ D);
       D'*(M \ n), D'*(M \ D)];

   % define linear objective vector
   c = [n'*(v + ha); D'*(v + ha)];

   % define inequality constraint matrix
   A = [n'*(M \ n), n'*(M \ D);
        mu -E'];

   % define inequality constraint vector
   b = [-n'*(v + ha); 0];

   % setup the LCP matrix and vector
   M_hat = [G -A'; A zeros(2)];
   q_hat = [c; -b]; 
   z = LCP(M_hat,q_hat);

   % now add in normal restitution
   z(1) = z(1) * (1 + e);

   % make sure that solution is still ok
   v_plus = [M \ n, M \ D]*[z(1:3)]+(v + ha);
   if (n'*v_plus < 0)
     % re-solve using new velocity
     c = [n'*v_plus; D'*v_plus];
     b = [-n'*v_plus; 0];
     q_hat = [c; -b]; 
     z = LCP(M_hat, q_hat);
   end
   
   %         if abs(z(1))<0.05
   %             e = 0;
   % %             vlplus1 = zeros(3,1)+1e-3;
   %         end
   


function fx = Example_3_Eqs_AD(Par)
% Equations of the benmark business cycle model with exogenous habits and a government spending shock and difference
% stationary growth. Equations are for algorithmic differention 
%
% Alfred Mauﬂner
% 
% 6 August 2019
%

import casadi.*;

% Parameters
astar=Par(1);
alpha=Par(2);
beta=Par(3);
eta=Par(4);
chi=Par(5);
delta=Par(6);
theta=Par(7);
psi=Par(8);
gstar=Par(9);

w=SX.sym('w',24);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
 k2=w(1);    k1=w(13);
ch2=w(2);   ch1=w(14);
 a2=w(3);    a1=w(15); %#ok<*NASGU>
 y2=w(4);    y1=w(16);
 c2=w(5);    c1=w(17);
 i2=w(6);    i1=w(18);
 N2=w(7);    N1=w(19);
 w2=w(8);    w1=w(20);
 r2=w(9);    r1=w(21);
 l2=w(10);   l1=w(22);
 z2=w(11);   z1=w(23);
 g2=w(12);   g1=w(24);

% equations of the model
f=[((c1-chi*ch1)^(-eta))*((1-N1)^(theta*(1-eta)))-l1;
   theta*(c1-chi*ch1)-(1-N1)*w1;
   w1-(1-alpha)*(y1/N1);
   r1-alpha*(y1/k1);
   y1-((a1*N1)^(1-alpha))*(k1^alpha);
   y1-c1-i1-exp(g1)*gstar;
   a1-exp(z1)*astar;
   a1*k2-(1-delta)*k1-i1;
   beta*(a1^(-eta))*l2*(1-delta+r2)-l1;
   ch2-(astar/a1)*((1-psi)*c1+psi*ch1);
   ];

fx=Function('HAB',{w},{f});
return;
end


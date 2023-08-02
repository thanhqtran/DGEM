function fx = Example_4_Eqs_AD(Par)
% Equations of the Calvo model.
%
% Alfred Mauﬂner
% 
% 16 August 2019
%
% The Calvo model is presented in Section 14.4.3 of my lecture note Computational Macroeconomics.

import casadi.*;

% Parameters
alpha=Par(1);
beta1=Par(2);
epsilon=Par(3);
eta=Par(4);
delta=Par(5);
delta1=Par(6);
delta2=Par(7);
gam0=Par(8);
gam1=Par(9);
phi=Par(10);
qstar=Par(11);
theta0=Par(12);
theta1=Par(13);
pistar=Par(14);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
s=SX.sym('w',44);

  k2=s(1);      k1=s(23);
  m2=s(2);	    m1=s(24);  
 qp2=s(3);     qp1=s(25);
pip2=s(4);    pip1=s(26); 
 sp2=s(5);     sp1=s(27);
  y2=s(6);	    y1=s(28);
  c2=s(7);		c1=s(29);
  i2=s(8);		i1=s(30);
  n2=s(9);		n1=s(31);
  w2=s(10);		w1=s(32);
  r2=s(11);		r1=s(33);
  l2=s(12);		l1=s(34);
 pi2=s(13);	   pi1=s(35);
 mu2=s(14);    mu1=s(36);
  q2=s(15);		q1=s(37);
 pa2=s(16);	   pa1=s(38);
  g2=s(17);     g1=s(39);
ga12=s(18);   ga11=s(40);
ga22=s(19);   ga21=s(41);
  s2=s(20);	    s1=s(42);
  z2=s(21);	    z1=s(43);
eps2=s(22);    eps1=s(44);

f=[(c1^(-eta))-l1;  
   theta0*(n1^theta1)-w1*l1;
   s1*y1-exp(z1)*(n1^(1-alpha))*(k1^alpha);
   w1-g1*(1-alpha)*exp(z1)*(n1^(-alpha))*(k1^alpha);
   r1-g1*alpha*exp(z1)*(n1^(1-alpha))*(k1^(alpha-1));
   y1-c1-i1;
   q1-(qstar^(1-delta1))*(qp1^delta1)*((pi1/pistar)^delta2)*exp(eps1);
   (((mu1/pi1)*m1)^(-gam1))-(l1/gam0)*((q1-1)/q1);  
   pa1-((epsilon/(epsilon-1))*ga11)/(pi1*ga21);
   1-(1-phi)*(pa1^(1-epsilon)) - phi*(pip1/pi1)^(1-epsilon);  
   s1-(1-phi)*(pa1^(-epsilon))-phi*((pip1/pi1)^(-epsilon))*sp1;   
   k2-(1-delta)*k1 - i1;
   m2-(mu1/pi1)*m1;
   sp2-s1;
   pip2-pi1;
   beta1*l2*(1-delta+r2)-l1;
   beta1*(q1/pi2)*l2-l1;
   ga11-g1*l1*(pi1^epsilon)*y1-beta1*phi*ga12;
   ga21-l1*(pi1^(epsilon-1))*y1-beta1*phi*ga22;
   qp2-q1;
   ];

   fx=Function('Calvo',{s},{f});
return;
end


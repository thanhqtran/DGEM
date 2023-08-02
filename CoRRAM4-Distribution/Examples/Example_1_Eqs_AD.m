function fx = Example_1_Eqs_AD(Par)
% Equations of the model of Example_1_AD.m for algorithmic derivatives
%
% Alfred Mauﬂner
% 
% 6 August July 2019
%

import casadi.*;

% Parameters
a=Par(1);
alpha=Par(2);
beta=Par(3);
eta=Par(4);
delta=Par(5);
theta=Par(6);

v=SX.sym('v',12);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
k2=v(1);    k1=v(7);
y2=v(2);    y1=v(8);
c2=v(3);    c1=v(9);
i2=v(4);    i1=v(10);
n2=v(5);    n1=v(11);
z2=v(6);    z1=v(12); %#ok<NASGU>

% equations of the model
f=[y1-exp(z1)*(n1^(1-alpha))*(k1^alpha);
   theta*c1-(1-n1)*(1-alpha)*(y1/n1);
   y1-c1-i1;
   a*k2-(1-delta)*k1-i1;
   1-beta*(a^(-eta))*((c1/c2)^eta)*(((1-n2)/(1-n1))^(theta*(1-eta)))*...
    (1-delta+alpha*(y2/k2));
    ];
 fx=Function('BM',{v},{f});

return;
end


function [fx, s] = Example_1_Eqs_SD()
% Equations of the model of Example_1_AD.m for symbolic derivatives
%
% Alfred Mauﬂner
% 
% 6 August July 2019
%

syms a eta delta theta;
beta=sym('beta');
alpha=sym('alpha');

syms k2 y2 c2 i2 n2 z2 k1 y1 c1 i1 n1 z1;

s=[k2 y2 c2 i2 n2 z2 k1 y1 c1 i1 n1 z1];

% equations of the model
fx=[y1-exp(z1)*(n1^(1-alpha))*(k1^alpha);
   theta*c1-(1-n1)*(1-alpha)*(y1/n1);
   y1-c1-i1;
   a*k2-(1-delta)*k1-i1;
   1-beta*(a^(-eta))*((c1/c2)^eta)*(((1-n2)/(1-n1))^(theta*(1-eta)))*...
    (1-delta+alpha*(y2/k2));
    ];

return;
end


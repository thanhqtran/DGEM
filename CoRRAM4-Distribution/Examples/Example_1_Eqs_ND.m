function f = Example_2_Eqs_ND(v,eqno,Par)
% Equations of the model of Example_1_ND.m for numeric derivatives
%
% Alfred Mauﬂner
% 
% 2 August July 2019
%

% Parameters
a=Par(1);
alpha=Par(2);
beta=Par(3);
eta=Par(4);
delta=Par(5);
theta=Par(6);

f=ones(5,1);

% variables of the model, 1 refers to period t and 2 to period t+1 variables
k2=v(1);    k1=v(7);
y2=v(2);    y1=v(8);
c2=v(3);    c1=v(9);
i2=v(4);    i1=v(10);
n2=v(5);    n1=v(11);
z2=v(6);    z1=v(12); %#ok<NASGU>

% equations of the model
f(1)=y1-exp(z1)*(n1^(1-alpha))*(k1^alpha);
f(2)=theta*c1-(1-n1)*(1-alpha)*(y1/n1);
f(3)=y1-c1-i1;
f(4)=a*k2-(1-delta)*k1-i1;
f(5)=1-beta*(a^(-eta))*((c1/c2)^eta)*(((1-n2)/(1-n1))^(theta*(1-eta)))*...
    (1-delta+alpha*(y2/k2));
if eqno~=0; f=f(eqno); end

return;
end


function f = Example_2_Eqs_ND(s,eqno,Par)
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
k2=s(1);    k1=s(8);
a2=s(2);    a1=s(9);
y2=s(3);    y1=s(10);
c2=s(4);    c1=s(11);
i2=s(5);    i1=s(12);
n2=s(6);    n1=s(13);
z2=s(7);    z1=s(14); %#ok<NASGU>

% equations of the model
f(1)=y1-((a*exp(z1)*n1)^(1-alpha))*(k1^alpha);
f(2)=theta*c1-(1-n1)*(1-alpha)*(y1/n1);
f(3)=y1-c1-i1;
f(4)=a1-a*exp(z1);
f(5)=a1*k2-(1-delta)*k1-i1;
f(6)=1-beta*(a1^(-eta))*((c1/c2)^eta)*(((1-n2)/(1-n1))^(theta*(1-eta)))*...
    (1-delta+alpha*(y2/k2));
if eqno~=0; f=f(eqno); end

return;
end


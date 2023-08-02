function [cvec,lvec,nvec,kvec,ParOut] = OLG_BGP1(x,Par)
%Returns the age profiles of consumption, marginal utility, labor supply, and wealth
%{
    Alfred Maußner

    2 September 2019, first version

    The function computes the age profiles of consumption, marginal utility,
    labor supply, and capital of a T-generation OLG model for the set of
    parameters given in Par and the aggregate stock of capital, labor, and
    the capital stock of the T-year old agents in x.

    In addition to those profiles it returns serveral aggregate results in ParOut.

%}

alpha=Par.alpha;
beta=Par.beta;
eta=Par.eta;
delta=Par.delta;
theta=Par.theta;
zeta=Par.zeta;
T=Par.T;
R=Par.R;

K=x(1);
L=x(2);
kT=x(3);

% arrays
cvec=zeros(T,1);
kvec=zeros(T,1);
nvec=zeros(R-1,1);
lvec=zeros(T,1);

% compute w and r
w=(1-alpha)*(K^alpha)*(L^(-alpha));
r=alpha*(K^(alpha-1))*(L^(1-alpha));

% compute b
q=(T-(R-1))/(R-1);
tau=(q*zeta)/(1+q*zeta);
b=(1-tau)*w*zeta*(T/(R-1))*L;

% compute sequence of Lagrange multipliers
zf=1-delta+r;

bR=beta*zf;
cvec(T)=zf*kT+b;
lvec(T)=cvec(T)^(-eta);

for s=T:-1:2
    lvec(s-1)=bR*lvec(s);
end

for s=(T-1):-1:R
    cvec(s)=lvec(s)^(-1/eta);
end

tmp1=theta/((1-tau)*w);
tmp2=tmp1^((theta*(eta-1))/(theta*(1-eta)-eta));

for s=R-1:-1:1
    cvec(s)=tmp2*lvec(s)^(1/(theta*(1-eta)-eta));
    nvec(s)=1-tmp1*cvec(s);
end
%nvec(nvec>0)=0;

kvec(T)=kT;
for s=T-1:-1:R
    kvec(s)=(kvec(s+1)+cvec(s)-b)/zf;
end

for s=R-1:-1:1
    kvec(s)=(kvec(s+1)+cvec(s)-(1-tau)*w*nvec(s))/zf;
end

ParOut.C=sum(cvec)/T;
ParOut.L=sum(nvec)/T;
ParOut.K=sum(kvec)/T;
ParOut.Y=(K^alpha)*(L^(1-alpha));
ParOut.I=ParOut.Y-ParOut.C;
ParOut.b=b;
ParOut.tau=tau;
ParOut.w=w;
ParOut.r=r;

return;

end


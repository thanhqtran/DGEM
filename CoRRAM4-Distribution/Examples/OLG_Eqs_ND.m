function fx = OLG_Eqs_ND(svec,eqno,Par)
%System of equations of the T generation OLG model
%{
   Alfred Mauﬂner

   3 September 2019

   Returns the system of 3T+R+5 equations that defines the dynamics of the T
   generation OLG model.

%}

% parameters
alpha=Par.alpha;
beta=Par.beta;
delta=Par.delta;
eta=Par.eta;
theta=Par.theta;
T=Par.T;
R=Par.R;
tau=Par.tau;
b=Par.b;

% recover variables from svec
nx=T-1;
ny=2*T+R+6;
nz=1;
nvar=nx+ny+nz;

kvec2=svec(1:T-1);
kvec1=svec(nvar+1:nvar+T-1);

Y1=svec(nvar+T);
C1=svec(nvar+T+1);
I1=svec(nvar+T+2);
K1=svec(nvar+T+3);
L1=svec(nvar+T+4);
w1=svec(nvar+T+5);
r2=svec(T+6);
r1=svec(nvar+T+6);
cvec1=svec(nvar+T+7:nvar+2*T+6);
nvec1=svec(nvar+2*T+7:nvar+2*T+R+5);
lvec2=svec(2*T+R+6:3*T+R+5);
lvec1=svec(nvar+2*T+R+6:nvar+3*T+R+5);
z1=svec(2*nvar);

fx=ones(nx+ny,1);

fx(1)=K1-sum(kvec1)/T;
fx(2)=L1-sum(nvec1)/T;
fx(3)=C1-sum(cvec1)/T;
fx(4)=Y1-exp(z1)*(K1^alpha)*(L1^(1-alpha));
fx(5)=w1-(1-alpha)*(Y1/L1);
fx(6)=r1-alpha*(Y1/K1);
fx(7)=Y1-C1-I1;
fx(8:R+6)=lvec1(1:R-1)-(cvec1(1:R-1).^(-eta)).*((1-nvec1).^(theta*(1-eta)));
fx(R+7:T+7)=lvec1(R:T)-(cvec1(R:T).^(-eta));
fx(T+8:T+R+6)=theta*cvec1(1:R-1)-(1-tau)*w1*(1-nvec1);
fx(T+R+7)=kvec1(T-1)*(1-delta+r1)+b-cvec1(T);
fx(T+R+8)=kvec2(1)+cvec1(1)-(1-tau)*w1*nvec1(1);
fx(T+R+9:T+2*R+6)=kvec2(2:R-1)+cvec1(2:R-1)-(1-delta+r1)*kvec1(1:R-2)-(1-tau)*w1*nvec1(2:R-1);
fx(T+2*R+7:2*T+R+6)=kvec2(R:T-1)+cvec1(R:T-1)-b-(1-delta+r1)*kvec1(R-1:T-2);
fx(2*T+R+7:3*T+R+5)=lvec1(1:T-1)-beta*lvec2(2:T)*(1-delta+r2);

if eqno~=0; fx=fx(eqno); end
return;
end


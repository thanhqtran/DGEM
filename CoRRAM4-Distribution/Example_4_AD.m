% Example_4_AD.m
%
%{
 
    Author:         Alfred Mauﬂner
    First Version:  12 March 2019, from Calvo_1.m
                    16 August 2019, extended to 3rd order solution via automatic differentiation

    Purpose:

    Example script to solve the Calvo model (as presented in my lecture notes,
    Computational Macroeconomics, in Section  14.4.3. The code to implement this
    model is taken from my Gauss program Money_5b.g. The equations of this model
    are set up in the function Calvo_1_Eqs_ND.m for numeric differentiation.

%}
clearvars;


% Parameters
alpha=0.27;
beta1=0.994;
eta=2.0;
delta=0.011;
epsilon=6;
nstar=0.13;
mustar=1.0171;
theta1=2.0;
epsmi=-0.2;
cm=0.837;
phi=0.25;
rhoz=0.90;
sigmaz=0.0072;
rhoeps=0.00;
sigmaeps=0.01;
delta1=0.8; 
delta2=1.5; 

% Stationary solution
gstar=(epsilon-1)/epsilon;
pastar=1;
sstar=1;
yk=(1-beta1*(1-delta))/(beta1*alpha*gstar);
kn=yk^(1/(alpha-1));
kstar=kn*nstar;
ystar=yk*kstar;
istar=delta*kstar;
cstar=ystar-istar;

pistar=mustar;
qstar=pistar/beta1;

wstar=gstar*(1-alpha)*(ystar/nstar);
rstar=gstar*alpha*(ystar/kstar);
lstar=cstar^(-eta);

theta0=wstar*lstar*(nstar^(-theta1));

mstar=mustar*(cstar/cm);
gam1=-(1/epsmi)*(1/qstar);

gam0=((qstar-1)/qstar)*(cstar^(-eta))*(mstar^gam1);

gamma1star=(gstar*lstar*(pistar^epsilon)*ystar)/(1-beta1*phi);
gamma2star=(lstar*(pistar^(epsilon-1))*ystar)/(1-beta1*phi);

s0=[kstar;mstar;qstar;pistar;sstar;ystar;cstar;istar;nstar;wstar;rstar;lstar;pistar;mustar;qstar;pastar;gstar; ...
   gamma1star;gamma2star;sstar];
% the same ordering as in Calvo_1.m for CoRRAM-M version 1
%x1=[kstar;mstar;qstar;pistar;sstar;0;0;ystar;cstar;istar;nstar;wstar;rstar;sstar;pistar;mustar;qstar;pastar;gstar; ...
%   gamma1star;gamma2star;lstar];

% create model environment
%Symbols={'k','m','qp','pip','sp', 'z', 'eps', 'y', 'c', 'i','n', 'w', 'r', 'l', 'pi', 'mu', 'q', ...
%         'pa', 'g', 'ga1', 'ga2', 's'}; 

%Symbols1={'k','m','qp','pip','sp', 'z', 'eps', 'y', 'c', 'i','n', 'w', 'r', 's', 'pi', 'mu', 'q', ...
%         'pa', 'g', 'ga1', 'ga2', 'l'}; 

%Calvo=DSGE(5,15,2,11,x0,Symbols);
%Calvo=DSGE(5,15,2,11,x0); % in case of numeric differentiation
Calvo=DSGE(5,15,2,11,s0);

Calvo.ParVal=[alpha;beta1;epsilon;eta;delta;delta1;delta2;gam0;gam1;phi;qstar;theta0;theta1;pistar];
Calvo.ParSym={'alpha','beta1','epsilon','eta','delta','delta1','delta2','gam0','gam1',...
              'phi','qstar','theta0','theta1','pistar'};
%Calvo.Equations='Calvo_1_Eqs_ND';
Calvo.Derivatives='AD';
Calvo.Equations='Example_4_Eqs_AD';

Calvo.Rho(1,1)=rhoz;
Calvo.Rho(2,2)=rhoeps;
Calvo.Omega(1,1)=sigmaz;
Calvo.Omega(2,2)=sigmaeps;

Calvo.Outfile='Example_4_AD';
Calvo.order=1;
Calvo.Flags.Reduce=0;
Calvo.Flags.Balance=1;
%Calvo.Flags.LoadDHT=1;

% solve Model
%Calvo.Skew=[[sigmaz 0 0 0];[0 0 0 sigmaeps]];
Calvo=SolveModel(Calvo);
if Calvo.rc~=0;
    error(Calvo.Messages{Calvo.rc});
end
error('stoped by user');

% simulation
analytic=1;
if ~analytic
    Calvo.VarBnd(:,1)=0.50*Calvo.VarVal;
    Calvo.VarBnd(:,2)=1.50*Calvo.VarVal;
    Calvo.Flags.CheckBounds=1;
    Calvo.Flags.Simulate=1;
    Calvo.nobs=120;
    Calvo.nofs=500;
    Calvo.Flags.Confidence=1;    
else
    Calvo.Filter='None';
    Calvo.Flags.Simulate=0;
    Calvo.Flags.Sylvester=0;
end

Calvo=SimulateModel(Calvo);    

% information for graphs and tables
Namestr={'Output';'Consumption';'Investment';'Hours';'Inflation';'TFP shock';'Interest rate shock'};
plotinfo=[[6 2];[7 2];[8 2];[9 3];[13 4];[21 1];[22 1]];
printidx=[6;7;8;9;13];
corridx=[6;13];

MakeGraphs(Calvo,'Charter',Namestr,plotinfo);

MakeTable(Calvo,Namestr(1:5),printidx,corridx);

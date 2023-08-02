% Example_3_AD.m
%
%{
 
    Author:         Alfred Mauﬂner
    First Version:  6 August 2019 (adapted from Habits_2.m)

    Purpose:

    Example script to solve the benchmark buiness cycle model with an exogenous
    consumption habit and government spending shocks.

    
%}

clearvars;

% Parameters of the model
a=exp(0.0022);
alpha=0.37;
beta=0.994;
eta=2.0;
chi=0.4;
delta=0.011;
Nstar=0.123;
psi=0.45;
Rho_z=0;
Sigma_z=0.0076;
gam=0.21;
Rho_g=0.55;
Sigma_g=0.0094;


% stationary solution
yk=(a^eta-beta*(1-delta))/(alpha*beta);
kstar=a*(yk^(1/(alpha-1)))*Nstar;
istar=(a-1+delta)*kstar;
ystar=yk*kstar;
gstar=gam*ystar;
cstar=ystar-istar-gstar;
chstar=cstar;
wstar=(1-alpha)*(ystar/Nstar);
rstar=alpha*yk;
theta=(1-Nstar)*(wstar/(cstar-chi*chstar));
lstar=((cstar-chi*chstar)^(-eta))*((1-Nstar)^(theta*(1-eta)));

% create model environment
HAB=DSGE(2,8,2,7);
HAB.VarVal=[kstar;chstar;a;ystar;cstar;istar;Nstar;wstar;rstar;lstar;0;0];
HAB.ParVal=[a;alpha;beta;eta;chi;delta;theta;psi;gstar];
HAB.VarXi =[1;1;0;0;0;1;1;1;0;1;0;-eta];
HAB.Equations='Example_3_Eqs_AD';
HAB.Rho(1,1)=Rho_z;
HAB.Rho(2,2)=Rho_g;
HAB.Omega(1,1)=Sigma_z;
HAB.Omega(2,2)=Sigma_g;
HAB.Outfile='Example_3_AD';
HAB.Derivatives='AD';
HAB.Flags.Sylvester=1;
HAB.Flags.Ds=1;
HAB.order=2;
HAB.inobs=15;
HAB.Loadpath=cd;

% solve Model
HAB=SolveModel(HAB);
if HAB.rc~=0
    error(HAB.Messages{HAB.rc});
end

HAB=SimulateModel(HAB);

Names={'Output','Consumption','Hours','TFP-Shock','Government spending shock'};
ipl=[[4 2];[5 2];[7 3];[11 2];[12 1]];
MakeGraphs(HAB,'Charter',Names,ipl);

HAB.Var(4).Name='Output';
HAB.Var(4).Plotno=2;

HAB.Var(5).Name='Consumption';
HAB.Var(5).Plotno=2;

HAB.Var(7).Name='Hours';
HAB.Var(7).Plotno=2;

HAB.Var(11).Name='TFP-Shock';
HAB.Var(11).Plotno=1;

HAB.Var(12).Name='Government Spending Shock';
HAB.Var(12).Plotno=1;

MakeGraphs(HAB,'Charter');


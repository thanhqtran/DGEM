%{

    Alfred Maußner

    2 August 2019

    Example_2_ND.m

    Purpose: Solve and simulate the model of Section 3.4.2 of the CoRRAM user manual.

%}

clearvars;

% Parameters of the model: preferences
beta =0.99; % discount factor
eta  =2.0;  % CRRA
theta=2.0; % utility weight of leisure

% production
a=1.004;    % growth factor
alpha=0.36; % capital share

% capital accumulation
delta=0.025; % rate of capital depreciation

% TFP shock
rho_nu=0;
sigma_nu=0.018;

% compute stationary solution
yk=(a^eta-beta*(1-delta))/(alpha*beta);
kn=a*yk^(1/(alpha-1));
ck=yk-(a-1+delta);
nstar=(theta/(1-alpha))*(ck/yk);
nstar=1/(1+nstar);
kstar=kn*nstar;
ystar=yk*kstar;
istar=(a-1+delta)*kstar;
cstar=ystar-istar;

% create an instance of the DSGE structure
nx=1;
ny=5;
nz=1;
nu=4;
BM=DSGE(nx,ny,nz,nu);

% supply information to the model
BM.Equations='Example_2_Eqs_ND';
BM.VarVal=[kstar;a;ystar;cstar;istar;nstar;0];
BM.ParVal=[a;alpha;beta;eta;delta;theta];
BM.Rho(1,1)=rho_nu;
BM.Omega(1,1)=sigma_nu;

BM.Derivatives='ND';
BM.order=2;

BM.Outfile='Example_2_ND';

BM=SolveModel(BM);
if BM.rc>0
    error(BM.Messages{BM.rc})
end

BM.Flags.Ds=1;
BM.VarXi=[1;0;1;1;1;0];
BM.burnin=0;
BM.Flags.Simulate=1;
BM.Flags.Confidence=1;
BM=SimulateModel(BM);

BM.Var(2).Name='Growth Factor';
BM.Var(2).Plotno=2;

BM.Var(3).Name='Output';
BM.Var(3).Print=1;
BM.Var(3).Corr=1;
BM.Var(3).Rel=1;
BM.Var(3).Plotno=3;

BM.Var(4).Name='Consumption';
BM.Var(4).Print=1;
BM.Var(4).Plotno=3;

BM.Var(5).Name='Investment';
BM.Var(5).Print=1;
BM.Var(5).Plotno=3;

BM.Var(6).Name='Hours';
BM.Var(6).Print=1;
BM.Var(6).Plotno=4;

BM.Var(7).Name='TFP Shock';
BM.Var(7).Plotno=1;

%Graphs(BM,'Charter');

%Tables(BM);
BM.Flags.Excel=1;
%MakeTable(BM);

Names={'Output','Consumption','Investment','Hours'};
ip=[3;4;5;6];
ic=[3;6];
MakeTable(BM,Names,ip,ic);


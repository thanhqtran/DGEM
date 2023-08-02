%{

    Alfred Maußner

    2 August 2019

    Example_1_ND.m

    Purpose: Solve and simulate the model of Section 3.4.1 of the CoRRAM user manual. The algorithm employs numeric differentiation.

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
rhoz=0.95;
sigmaz=0.007;

% compute stationary solution
yk=(a^eta-beta*(1-delta))/(alpha*beta);
kn=yk^(1/(alpha-1));
ck=yk-(a-1+delta);
nstar=(theta/(1-alpha))*(ck/yk);
nstar=1/(1+nstar);
kstar=kn*nstar;
ystar=yk*kstar;
istar=(a-1+delta)*kstar;
cstar=ystar-istar;

% create an instance of the DSGE structure
nx=1;
ny=4;
nz=1;
nu=3;
BM=DSGE(nx,ny,nz,nu);

% supply information to the model
BM.Equations='Example_1_Eqs_ND';
BM.VarVal=[kstar;ystar;cstar;istar;nstar;0];
BM.ParVal=[a;alpha;beta;eta;delta;theta];
BM.Rho(1,1)=rhoz;
BM.Omega(1,1)=sigmaz;

BM.Derivatives='ND';
BM.order=2;

BM.Outfile='Example_1_ND';

BM=SolveModel(BM);
if BM.rc>0
    error(BM.Messages{BM.rc})
end

BM=SimulateModel(BM);

% version one to produce graphs and tables
BM.Var(2).Name='Output';
BM.Var(2).Print=1;
BM.Var(2).Corr=1;
BM.Var(2).Plotno=2;

BM.Var(3).Name='Consumption';
BM.Var(3).Print=1;
BM.Var(3).Plotno=2;

BM.Var(4).Name='Investment';
BM.Var(4).Print=1;
BM.Var(4).Plotno=2;

BM.Var(5).Name='Hours';
BM.Var(5).Print=1;
BM.Var(5).Plotno=2;

BM.Var(6).Name='TFP-Shock';
BM.Var(6).Plotno=1;

MakeGraphs(BM,'Arial');
MakeTable(BM);

% version two to produce graphs and tables
Names={'Output','Consumption','Investment','Hours','TFP-Shock'};
ipl=[[2 2];[3 2];[4 2];[5 2];[6 1]];
MakeGraphs(BM,'Arial',Names,ipl);

ipr=[2;3;4;5];
ic=[2];
MakeTable(BM,Names(1:4),ipr,ic);

%{

    Alfred Maußner

    2 August 2019, first version
   25    May 2022, recent version

    Example_1_AD.m

    Purpose: Solve and simulate the model of Section 3.4.1 of the CoRRAM user manual. This version employs automatic differention.

    Remarks: 

        1) be shure to have installed the CasAD toolbox before you run this script and adjust the path in line 17
        2) different from the script Example_1_ND.m this script employs a simulation with the order of the solution
           to compute impulse response. (The MakeGraphs function always employs the first-order solution!) 

%}

clearvars;
addpath('d:\code\othercode\CasAD\casadi-windows-matlabR2014b-v3.4.5'); % path to the CasADi files, change according to your local directory structure

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
BM.Equations='Example_1_Eqs_AD';
BM.VarVal=[kstar;ystar;cstar;istar;nstar;0];
BM.ParVal=[a;alpha;beta;eta;delta;theta];
BM.Rho(1,1)=rhoz;
BM.Omega(1,1)=sigmaz;

BM.Derivatives='AD';
BM.order=2;

BM.Outfile='Example_1_AD';
BM.Loadpath=cd;
BM=SolveModel(BM);
if BM.rc>0
    error(BM.Messages{BM.rc})
end

% compute IRFs (added 26 May 2022)
BM.burnin=1000;
BM.inobs=16;
vhut = IRF_Sim(BM,1,[]);

% plot
Data(1).x=[vhut(:,2),vhut(:,6)]; Data(1).n={'Output','TFP Shock'};
Data(2).x=vhut(:,3); Data(2).n='Consumption';
Data(3).x=vhut(:,4); Data(3).n='Investment';
Data(4).x=vhut(:,5); Data(4).n='Hours';
nr=2; nc=2;
fname='Century Schoolbook';
fsize=12;
fmult=1.5;
Plot1(Data,nr,nc,fname,fsize,fmult);

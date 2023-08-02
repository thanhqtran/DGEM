%{
    OLG.m

    Author:         Alfred Maußner
    First Version:  2 March 2020

    Purpose:

    Example script to solve and simulate the overlapping generations model in
    Heer and Maußner, General Equilibrium Modeling, 2nd Edition, Springer: 2009
    with the Matlab version of my toolbox CoRRAM. The script employs numeric differentiation.

%
    
%}

clearvars;
addpath('d:\code\mycode\matlab\corram4'); % change the path to your installation of CoRRAM!
Par.alpha=0.3;
Par.beta=0.99;
Par.eta=2.0;
Par.delta=0.04;
Par.theta=2.0;
Par.zeta=0.3;
Par.T=240;
Par.R=161;

rho_z=0.814;
sigma_z=0.0142;

x0=[5;0.45;1];

fun=@(z)OLG_GetBGP(z,Par);
foptions=optimoptions('fsolve');
foptions.FiniteDifferenceType='central';
foptions.FunctionTolerance=1e-10;
[x1,fx,rc]=fsolve(fun,x0,foptions);

if rc<0
    error('Not able to compute balanced growth path');
end

%Plot 
[cvec,lvec,nvec,kvec,ParOut] = OLG_BGP1(x1,Par);
data(1).y=kvec; data(1).x=(1:Par.T)'; data(1).t='Capital stocks';
data(2).y=cvec; data(2).x=(1:Par.T)'; data(2).t='Consumption';
data(3).y=nvec; data(3).x=(1:Par.R-1)'; data(3).t='Labor supply';

myColor=[[0 0 0];[0 0 1];[0 0.44 0];[1 0 0];[0 0.49 0.49];[0.46 0 0.46];[0.58 0.29 0];[0 0 0]];
myLineStyle={'-','-','-','-','-','-','-','--','--','--','--','--'};

noplot=1;
if ~noplot
    for pn=1:3
        subplot(2,2,pn);
        hl=plot(data(pn).x,data(pn).y);
        ax=gca;
        ax.FontName='Charter';
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier=1.2;   
        xlabel('Generation');
        xu=max(data(pn).x);
        ax.XLim=[1 xu];
        hl(1).LineStyle=myLineStyle{1};
        hl(1).LineWidth=1;
        hl(1).Color=myColor(1,:);
        title(data(pn).t,'FontName','Charter');
    end
end
Par.tau=ParOut.tau;
Par.b  =ParOut.b;

%svec=[kvec(2:Par.T);ParOut.Y;ParOut.C;ParOut.I;ParOut.K;ParOut.L;ParOut.w;ParOut.r;cvec;nvec;lvec;0];

%fx=OLG_Eqs_ND([svec;svec],0,Par);

nx=Par.T-1;
ny=2*Par.T+Par.R+6;
nz=1;
nu=Par.T+Par.R+7;
svec=[kvec(2:Par.T);ParOut.Y;ParOut.C;ParOut.I;ParOut.K;ParOut.L;ParOut.w;ParOut.r;cvec;nvec;lvec];
OLG=DSGE(nx,ny,nz,nu,svec);
OLG.ParVal=Par;
OLG.Rho=rho_z;
OLG.Omega=sigma_z;
OLG.order=1;
OLG.etol=1e-8;
OLG.Derivatives='ND';
OLG.Equations='OLG_Eqs_ND';
OLG.Outfile='OLG_1';
OLG=SolveModel(OLG);
if OLG.rc>0
    error(OLG.Messages(OLG.rc));
end
OLG.Flags.Simulate=0; % analytical second moments
OLG.hpl=100;          % weight of HP-filter
OLG=SimulateModel(OLG);
nvar=nx+ny+nz;
Namestr={'Output';'Consumption';'Investment';'Hours';'Real Wage';'Log of TFP Shock'};
plotinfo=[[nx+1 2];[nx+2 2];[nx+3 2];[nx+5 3];[nx+6 3];[nvar 1]];
MakeGraphs(OLG,'Charter',Namestr,plotinfo);
printidx=[nx+1;nx+2;nx+3;nx+5;nx+6];
corridx=[nx+1;nx+5];
MakeTable(OLG,Namestr(1:5),printidx,corridx);
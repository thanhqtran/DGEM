function vhut = IRF_Sim(Mod,shock,PF)
% Computes impulse responses to a shock in a model with stochastic growth of labor augmenting technical progress
%{
 
  Alfred Maußner

  12 January 2021, IRF_SG_Proj.m
  26     May 2022, renamed to IRF_Sim.m and generalized to arbitray DSGE models

  Purpose: given the solution of an instance of the DSGE classs in Mod, compute the time path
           of the model's state variables x(t+1) and the jump variables y(t) from the policy functions PF
 
           [xt(:,t+1), yt(:,t)]=PF(v),v=[x(:,t);zt(:,t)])  

           where

           zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*eps(:,t+1), eps(:,t) iid N(0,eye(nz)).

           The policy function PF can be any admissible function that returns the same
           set of variables [xt(:,t+1) yt(:,t)] as the policy function Mod.PFXY! PF is the
           pointer to this Matlab function. If PF is empty the function employs the policy
           function from the perturbation solution defined in Mod.PF. 

           Since for non-linear solutions the solution usually does not return to the
           deterministic steady state, the function first computes the stochastic steady state by
           iterating for Mod.burnin periods over the above equation starting with xt(:,1)=xstar
           and zt(:,t)=0 for t=1, 2, ..., Mod.burin.

           The values which the variables attain after Mod.burin periods represent the stochastic stationary
           solution.

           In period tb+1 (where tb=Mod.burnin), there is a shock in the variable zt(shock,tb+1) of size Mod.Omega(shock,shock)
           and the iterations continue to Mod.burnin+Mod.inobs.

           The impulse responses are defined as the log-deviations of the levels of the variables on the
           shocked time path and the un-shocked path. 

           In models without per-capita growth or in models with deterministic trend in labor augmenting technical progress A(t)
           the log-deviations of the variables in levels are equal to the log-deviations of the stationary variables (either
           scaled variables as output or unscaled as working hours). In models with per capita growth (either because the
           growth factor is endogenous or follows an exogenous stochastic process) the shocked path of a variable i in levels is
           given by:
            
           V(i,t)=A(t-1)^xi(i)*v(i,t), A(t)=yt(1,t)*A(t-1), A(tb)=1, t=tb+tau

           and the unshocked path is given by

           VT(i,t)=(yt(1,tb)^(tau-1))*v(i,tb)

           for v(:,t):=[xt(:,t);yt(:,t)]; i=1, ..., Mod.nx+Mod.ny.

           Note that the function assumes that the first element in the vector of jump variables yt(:,t) is the growth factor a(t):=A(t)/A(t-1)

           If the flag Mod.Flags.DS is true, the function draws information on the scaling factors xi(i) from the vector
           Mod.Xi. Otherwise the function assumes a model without growth or deterministic growth factor a.

 
           The function returns deviations in percent!

   Reference: Alfred Maußner, Impulse Response Functions from Non-Linear Solutions of DSGE Models, University
   of Augsburg, 2021

%}

% check if PF is empty
if isempty(PF)
    PF=@(v)Mod.PFXY(v);
end

tb=Mod.burnin;
ti=Mod.inobs;
nobs=tb+ti;
nx=Mod.nx;
nz=Mod.nz;
ny=Mod.ny;

xt=zeros(nx,nobs+1);
zt=zeros(nz,nobs+1);
yt=zeros(ny,nobs);
xstar=Mod.VarVal(1:Mod.nx); % that is where iterations start from
xt(:,1)=xstar;

% burnin to obtain stochastic stationary solution
for t=1:tb
    v=[xt(:,t);zt(:,t)]';
    [xt(:,t+1), yt(:,t)]=PF(v);
end
vbar=[xt(:,tb);yt(:,tb)];

% shock occurs in t=tb+1
et=zeros(nz,ti+1);
et(shock,1)=1;

% continue simulations
for t=tb:nobs
    zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*et(:,t-(tb-1));
    v=[xt(:,t);zt(:,t)]';
    [xt(:,t+1), yt(:,t)]=PF(v);
end

% log deviations
vhut=[log(xt(:,tb:nobs))-log(vbar(1:nx));log(yt(:,tb:nobs))-log(vbar(1+nx:nx+ny))]';

if Mod.Flags.Ds
    tmp=[0;cumsum(log(yt(1,tb:nobs-1))-log(vbar(nx+1)))'];
    xivec=Mod.VarXi';
    vhut=vhut+xivec.*tmp;
end
vhut=100*[vhut, zt(shock,tb:nobs)'];

return;

end


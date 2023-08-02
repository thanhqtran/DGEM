function [Ex,Ey,Ez] = First_Moments( Mod)
% Computes the unconditional expectations of the variables of a DSGE model

%{
    Copyright: Alfred Maußner

    Purpose: Computes the unconditional expected means of the variables of
             a DSGE model. See equations (4.15) and (4.16) of Heer and Maußner,
             Dynamic General Equilibrium Modeling, 3rd. Edition.

    Revision history:

        17 June 2019 (from a previous version of Moments_TD.m in CoRRAM-M)

    Input:  Mod, an instance of the DSGE class

    Output: Ex, nx by one vector, the means of the state variables,
            Ey, ny by one vector, the means of the jump variables
            Ez, nz by one vector, the means of the shocks

 %}

% for simplicity
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nw=nx+nz;

% compute Gamma^w_0
if nx>0    
    Sigma_nu=[zeros(nx,nx+nz);[zeros(nz,nx), Mod.Omega*Mod.Omega']];
    Amat=[Mod.Hxw;[zeros(nz,nx), Mod.Rho]];
else
    Sigma_nu=Mod.Omega*Mod.Omega';
    Amat=Mod.Rho;
end
Gamma0=Lyapunov(Amat,Sigma_nu);

% compute means
Bmat=eye(nw)-Amat;
if nx>0
    Hww=[Mod.Hxww;zeros(nz*nw,nw)];
    bvec=zeros(nw,1);
    for i =1:nw
        bvec(i)=sum(diag(Hww((i-1)*nw+1:i*nw,:)*Gamma0));
    end
    if ~isempty(Mod.Mu)
        bvec=0.5*bvec+0.5*[Mod.Hxss;zeros(nz,1)]+[zeros(nx,1);Mod.Mu];
    else
        bvec=0.5*bvec+0.5*[Mod.Hxss;zeros(nz,1)];
    end
    Ew=linsolve(Bmat,bvec);
    Ex=Ew(1:nx);
    Ez=Ew(1+nx:nw);
else
    Ex=[];
    if ~isempty(Mod.Mu)
        Ez=linsolve(Bmat,Mod.Mu);
    else
        Ez=zeros(nz,1);
    end
end
bvec=zeros(ny,1);
for i=1:ny
    bvec(i)=sum(diag(Mod.Hyww((i-1)*nw+1:i*nw,:)*Gamma0))+Mod.Hyss(i);
end
if nx>0 
    Ey=Mod.Hyw*Ew+0.5*bvec;
else
    Ey=Mod.Hyw*Ez+0.5*bvec;
end

return;

end


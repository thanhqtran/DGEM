function [FN,DN,HN,TN,rc] = GetFDHT_AD( Mod )
% Compute the Jacobian and if desired the Hessian and the matrix of third-order derivatives of the canonical DSGE model
%{
    Copyright:  Alfred Maußner    

    Purpose:  This function uses the CasADi toolbox to
              compute numerical derivatives of the system of equations that
              defines an instance of the canonical DSGE model (see DSGE.m)
              via automatic (or algorithmic) differentiation.

    Revision history:          
           
    8 March   2019, first version 
   15 March   2019, computation of Hessian and third-order matrix added
   27 May     2019,  ordering changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))

    Input: Mod, an instance of the DSGE class

    Output: FN, system of equations evaluated at the stationary solution
            DN, numeric Jacobian matrix
            HN, numeric Hessian matrix
            TN, numeric matrix of third-order derivatives
            rc, return code that defines errors

%}

import casadi.*;

% initialize
HN=[];
TN=[];
rc=0;

nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nxy=nx+ny;    
na=nxy+nz;
indc=[1:nx,1+nxy:na,1+nx:nxy,1+na:na+nx,1+na+nxy:2*na,1+na+nx:na+nxy];

% define file handle (existence is checked by invoking the DSGE class)
Sys=str2func(Mod.Equations);

% CasADi variable
W=SX.sym('W',2*(Mod.nx+Mod.ny+Mod.nz));

% system of equations
F_AD=Sys(Mod.ParVal);

% numeric value in Matlab 
FN=full(F_AD([Mod.VarVal;Mod.VarVal]));

if max(abs(FN))>Mod.etol % check stationary solution
    DN=[];
    rc=2;
    return;
end

% Jacobian
DF_AD=Function('DF',{W},{jacobian(F_AD(W),W)});
DN=full(DF_AD([Mod.VarVal;Mod.VarVal])); % evaluate to Matlab matrix
DN=DN(:,indc); % reordering of columns

% Check for invalid Jacobian
if any(any(isnan(DN))) || any(any(isinf(DN)))
    rc=3;
    return;
end

if Mod.order==1; return; end   

% compute Hessian
HF_AD=Function('HF',{W},{jacobian(DF_AD(W),W)});
H=full(HF_AD([Mod.VarVal;Mod.VarVal])); % evaluate to Matlab matrix

% arrange and reorder Hessian
HN=zeros(nxy*2*na,2*na);
indr=0:nxy:(2*na-1)*nxy;

for r=1:nxy
    tmp=H(r+indr,:);
    HN(1+(r-1)*2*na:r*2*na,:)=tmp(indc,indc);
end
if Mod.order<3; return; end

% compute matrix of third-order derivatives
TF_AD=Function('TF',{W},{jacobian(HF_AD(W),W)});
TT=full(TF_AD([Mod.VarVal;Mod.VarVal])); % evaluate to Matlab matrix

% arrange 
m=2*na;
TN=zeros(nxy*m*m,m);
indr1=(0:nxy:(m^2-1)*nxy)'; % extract m*m by m matrices
indr2=[1:m*nx, ...
      1+m*nxy:m*na, ...
      1+m*nx:m*nxy, ...
      1+m*na:m*(na+nx), ...
      1+m*(na+nx+ny):m*m, ...
      1+m*(na+nx):m*(na+nx+ny)];
indr3=repmat(indc,m,1)+(0:m:m*(m-1))';   
indr3=reshape(indr3',m^2,1);

for r=1:nxy
    indr1a=r+indr1;
    tmp=TT(indr1a,:);
    tmp=tmp(indr2,:);
    TN(1+(r-1)*m^2:r*m^2,:)=tmp(indr3,indc);
end

return;

end


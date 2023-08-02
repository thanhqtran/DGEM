function [Hxww, Hxss, Hyww, Hyss, rc]=Quadratic(Mod)
% Obtain quadratic part of the perturbation solution of the canonical DSGE models

%{

 Copyright: Alfred Maußner

 Purpose: Compute the coefficients of the matrices Hxww, Hxss, Hyww, and Hyss
          that form the quadratic part of the perturbation solution of the
          canonical DSGE model. This model is given by

          (CM1): 0=E_t g(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))
          (CM2): 0=z(t+1)-Rho*z(t)+sigma*Omega*epsilon(t+1).  

          The quadratic part of the perturbation solution at the point
 
          x(t)=x, y(t)=y, and z(t)=0

          for the vector of state variables x(t+1) is given by

          0.5*(I(nx) otimes wbar(t)')*Hxww*wbar(t) + 0.5*Hxss*sigma^2, where wbar(t):=[(x(t)-x);z(t)].

          For the vector of jump variables y(t) it is given by
 
          0.5*(I(ny) otimes wbar(t)')*Hyww*wbar(t) + 0.5*Hyss*sigma^2

          
 Input: Mod, an instance of the DSGE class that stores all matrices required to
        solve for the quadratic part.

 Ouput: Hxww, nx*nw by nw matrix, nx is the number of state variable, nz is the number of shocks, 
                          and nw=nx+nz

        Hxss, nx by 1 vector
        Hyww, ny*nw by nw matrix, ny is the number of jump variables
        Hyss, ny by 1 vector.

 Revisions history:
     06 December 2016, first version
     12 January  2017 (adapted to the DSGE class)
     11 May      2017 (non zero means of innovations)
     14 March    2019 from script to function
      8 November 2019 return values in case of nx==0 for Hxww and Hxss added
     24 November 2020 return values in case of nz==0 for Hxss and Hyss added

%}

% check whether first-order solution has been computed
if Mod.nx>0
    if isempty(Mod.Hyw) || isempty(Mod.Hxw); error('You must first solve for the linear part.'); end
else
    if isempty(Mod.Hyw); error('You must first solve for the linear part.'); end
end

if isempty(Mod.DN) || isempty(Mod.HN); error('The Jacobian and/or the Hessian matrices are missing.'); end

% initialize, for simplicty
nx=Mod.nx;
nz=Mod.nz;
ny=Mod.ny;
na=nx+ny+nz;
rc=0;
Hxss=[];
Hyss=[];

if nz>0
    if nx>0
        fw=[Mod.Hxw;...
           [zeros(nz,nx) Mod.Rho];...
           [Mod.Hyw(:,1:nx)*Mod.Hxw(:,1:nx) (Mod.Hyw(:,1:nx)*Mod.Hxw(:,nx+1:nx+nz)+Mod.Hyw(:,nx+1:nx+nz)*Mod.Rho)];...
           [eye(nx) zeros(nx,nz)];...
           [zeros(nz,nx) eye(nz)];
            Mod.Hyw];
    else
        fw=[Mod.Rho; Mod.Hyw*Mod.Rho;eye(nz);Mod.Hyw];
    end
else % the case nz==0 and nx== cannot occur
    fw=[Mod.Hxw;Mod.Hyw*Mod.Hxw;eye(nx);Mod.Hyw];
end
        
A1=kron(eye(nx+ny),fw')*Mod.HN*fw;
B1=kron(Mod.DN(:,1+2*na-ny:2*na),eye(nx+nz)); % gy
if nx>0; B2=kron(Mod.DN(:,1:nx),eye(nx+nz)); end % gx'
B3=kron(Mod.DN(:,1+nx+nz:nx+nz+ny),eye(nx+nz)); % gy'
if nz>0
    if nx>0
        C2=[Mod.Hxw;[zeros(nz,nx) Mod.Rho]];
    else
        C2=Mod.Rho;
    end
else
    C2=Mod.Hxw;
end
C1=kron(eye(ny),C2');
if nx>0; C3=kron(Mod.Hyw(:,1:nx),eye(nx+nz)); end
if Mod.Flags.Sylvester
    if nx>0
        atilde=[B1, (B2+B3*C3)];
    else
        atilde=B1;
    end
    btilde=-C2;
    ctilde=-A1;
    if nx>0
        dtilde=[B3*C1,zeros((nx+ny)*(nx+nz),nx*(nx+nz))];
    else
        dtilde=B3*C1;
    end
    etilde=eye(nx+nz);
    ftilde=zeros((nx+ny)*(nx+nz),nx+nz);
    rmat=Sylvester(atilde,btilde,ctilde,dtilde,etilde,ftilde);
    Hyww=rmat(1:ny*(nx+nz),:);
    if nx>0; Hxww=rmat(ny*(nx+nz)+1:(nx+ny)*(nx+nz),:); end
else
    bigmat=kron(eye(nx+nz),B1)+kron(C2',B3*C1);
    if nx>0; bigmat=[bigmat,kron(eye(nx+nz),(B2+(B3*C3)))]; end
    vecEG=linsolve(bigmat,-A1(:));
    Hyww=reshape(vecEG(1:ny*(nx+nz)^2),[ny*(nx+nz),nx+nz]);
    if nx>0; Hxww=reshape(vecEG(ny*(nx+nz)^2+1:(ny+nx)*(nx+nz)^2),[nx*(nx+nz),nx+nz]); end
end

if nz<0.5; Hxss=zeros(nx,1); Hyss=zeros(ny,1); return; end

% Solve for Hxss and Hyss
smat=Mod.Omega*(Mod.Omega');
if nx>0
    nmat=[zeros(nx,nz);eye(nz);Mod.Hyw(:,1+nx:nx+nz);zeros(nx+nz+ny,nz)];
else
    nmat=[eye(nz);Mod.Hyw(:,1+nx:nx+nz);zeros(nz+ny,nz)];
end
xvec=tracem(kron(eye(nx+ny),nmat')*Mod.HN*(nmat*smat));
M11=zeros(ny,1);
for ii=1:ny
   M11(ii,1)=trace(smat*Hyww((ii-1)*(nx+nz)+nx+1:ii*(nx+nz),1+nx:nx+nz));
end
xvec=xvec+(Mod.DN(:,1+nx+nz:na)*M11); 
if ~isempty(Mod.Muss) % this captures the case of shocks with the mean preserving spread property
    xvec=xvec+(Mod.DN(:,1+nx:nx+nz)+Mod.DN(:,1+nx+nz:na)*Mod.Hyw(:,1+nx:nx+nz))*Mod.Muss;
end
bigmat=Mod.DN(:,nx+nz+1:nx+nz+ny)+Mod.DN(:,2*(nx+nz)+ny+1:2*na);
if nx>0; bigmat=[bigmat,Mod.DN(:,1:nx)+Mod.DN(:,1+nx+nz:nx+nz+ny)*Mod.Hyw(:,1:nx)]; end
VecEG=linsolve(bigmat,-xvec);
Hyss=VecEG(1:ny);
if nx>0
    Hxss=VecEG(ny+1:ny+nx);
else
    Hxww=[];
    Hxss=[];
end

return;

end


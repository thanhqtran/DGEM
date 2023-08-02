function [Hxwww,Hywww,Hxssw,Hyssw,Hxsss,Hysss,rc] = Cubic(Mod)
%Computes the cubic part of the perturbation solution of the canonical DSGE model
%{
    Author: Alfred Maußner

    Purpose: Computes the coefficients of the matrices Hxwww, Hxssw, Hxsss,
             Hywww, Hyssw, and Hysss that form the cubic part of the 
             perturbation solution of the canonical DSGE model. This model is given by

          (CM1): 0=E_t g(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))
          (CM2): 0=z(t+1)-Rho*z(t)+sigma*Omega*epsilon(t+1).  

          The cubic part of the perturbation solution at the point
 
          x(t)=x, y(t)=y, and z(t)=0

          for the vector of state variables x(t+1) is given by

          (1/6)*(I(nx) otimes wbar(t)' otimes wbar(t)')*Hxwww*wbar(t) 
          + 0.5*(I(nx) otimes wbar(t)')*Hxssw*sigma^2 + (1/6)*Hxsss*sigma^3, 

          where wbar(t):=[(x(t)-x);z(t)].

          For the vector of jump variables y(t) it is given by
 
          (1/6)*(I(ny) otimes wbar(t)' otimes wbar(t)')*Hywww*wbar(t) 
         + 0.5*(I(ny) otimes wbar()')*Hyssw*sigma^2 + (1/6)*Hysss*sigma^3

    Input: Mod, an instance of the DSGE class that stores all matrices required to
          solve for the cubic part.

    Ouput: Hxwww, nx*nw^2 by nw matrix, nx is the number of state variable, nz is the number of shocks, 
                          and nw=nx+nz

           Hxssw, nx by nw matrix,
           Hxsss, nx by 1 vector,
           Hywww, ny*nw^2 by nw matrix
           Hyssw, ny by nw matrix
           Hysss, ny by 1 vector

    Revision history:
       
       27 March 2019, first version (from a previous version of the script Cubic.m)
       18 April 2019, bug, was Mod.Hywww instead of Hywww


%}

% initialize
Hxssw=[];
Hyssw=[];
Hxsss=[];
Hysss=[];
rc=0;

% check for required matrices
if isempty(Mod.Hxww) || isempty(Mod.Hyww); error('You must first solve for the quadratic part.'); end
if isempty(Mod.TN); error('The matrix of third-order derivatives has not been computed.'); end

% for simplicity: get required dimensions
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nw=nx+nz;
na=nx+ny+nz;

% Solve for Hxwww and Hywww
if nz>0
    fw=[Mod.Hxw;
        [zeros(nz,nx), Mod.Rho];
        [Mod.Hyw(:,1:nx)*Mod.Hxw(:,1:nx), Mod.Hyw(:,1:nx)*Mod.Hxw(:,1+nx:nx+nz)+Mod.Hyw(:,1+nx:nx+nz)*Mod.Rho];
        [eye(nx),zeros(nx,nz)];
        [zeros(nz,nx),eye(nz)];
        Mod.Hyw];
else
    fw=[Mod.Hxw;Mod.Hyw*Mod.Hxw;eye(nx);Mod.Hyw];
end
if nz>0
    uw=[Mod.Hxw;[zeros(nz,nx), Mod.Rho]];
else
    uw=Mod.Hxw;
end

if nz>0
    fww=[Mod.Hxww;
         zeros(nz*nw,nw);
         kron(eye(ny),uw')*Mod.Hyww*uw+kron(Mod.Hyw(:,1:nx),eye(nw))*Mod.Hxww;
         zeros(nx*nw,nw);
         zeros(nz*nw,nw);
         Mod.Hyww];
else
    fww=[Mod.Hxww;
         kron(eye(ny),uw')*Mod.Hyww*uw+kron(Mod.Hyw,eye(nw))*Mod.Hxww;
         zeros(nx*nw,nw);
         Mod.Hyww];
end

if nz>0
    uww=[Mod.Hxww;zeros(nz*nw,nw)];
else
    uww=Mod.Hxww;
end
Q1=kron(Mod.DN(:,1:nx),eye(nw^2));
Q2=kron(Mod.DN(:,1+nx+nz:na),eye(nw^2));
Q3=kron(Mod.DN(:,1+na+nx+nz:2*na),eye(nw^2));

P1=kron(kron(eye(nx+ny),fw'),fw')*Mod.TN*fw;
P2=kron(eye(nx+ny),Getfxxt(fww)')*Mod.HN*fw;
P3=kron(kron(eye(nx+ny),fw'),eye(nw))*kron(Mod.HN,eye(nw))*fww;
P4=kron(eye(nx+ny),Getfxt(fw)')*kron(Mod.HN,eye(nw))*fww;
P5=kron(eye(ny),Getfxxt(uww)')*Mod.Hyww*uw + ...
   kron(kron(eye(ny),uw'),eye(nw))*kron(Mod.Hyww,eye(nw))*uww + ...
   kron(eye(ny),Getfxt(uw)')*kron(Mod.Hyww,eye(nw))*uww;
P5=Q2*P5;

S1=kron(kron(eye(ny),uw'),uw');
S3=kron(Mod.Hyw(:,1:nx),eye(nw^2));

if Mod.Flags.Sylvester
   atilde=[Q3, (Q1+Q2*S3)];
   btilde=-uw;
   ctilde=-(P1+P2+P3+P4+P5);
   dtilde=[Q2*S1,zeros((nx+ny)*nw^2,nx*nw^2)];
   etilde=eye(nw);
   ftilde=zeros((nx+ny)*nw^2,nw);
   rmat=Sylvester(atilde,btilde,ctilde,dtilde,etilde,ftilde);
   Hywww=rmat(1:ny*nw^2,:);
   Hxwww=rmat(ny*nw^2+1:(nx+ny)*nw^2,:); 
else
    xvec=-(P1+P2+P3+P4+P5);    
    bigmat=[kron(eye(nw),Q3)+kron(uw',Q2*S1), ...
            kron(eye(nw),(Q1+Q2*S3))];
    vecEG=linsolve(bigmat,xvec(:));
    Hywww=reshape(vecEG(1:ny*nw^3),ny*nw^2,nw);    
    Hxwww=reshape(vecEG(1+ny*nw^3:(nx+ny)*nw^3),nx*nw^2,nw);    
end
if nz<0.5; return; end

% Hxssw and Hyssw, I recompute elements from Quadratic.m
smat=Mod.Omega*(Mod.Omega');
nmat=[zeros(nx,nz);eye(nz);Mod.Hyw(:,1+nx:nx+nz);zeros(nw+ny,nz)];
M11=zeros(ny,1);
for ii=1:ny
   M11(ii,1)=trace(smat*Mod.Hyww((ii-1)*(nx+nz)+nx+1:ii*(nx+nz),1+nx:nx+nz));
end
fss=[Mod.Hxss; zeros(nz,1); (M11+Mod.Hyss+Mod.Hyw(:,1:nx)*Mod.Hxss);zeros(nw,1);Mod.Hyss];
A1=kron(eye(nx+ny),fw')*Mod.HN*fss;
Hywzt=zeros(ny,nw*nz);
for i=1:ny
    Hywzt(i,:)=reshape(Mod.Hyww((i-1)*nw + (nx+1:nw),:),1,nw*nz);
end
fsw=[zeros(nw,nw*nz);(Hywzt*kron(uw,eye(nz)));zeros(nw+ny,nw*nz)];
A2=2*tracem(kron(eye(nx+ny),fsw')*Mod.HN*nmat*smat);
A3=tracem(kron(kron(eye(nx+ny),fw'),nmat')*Mod.TN*nmat*smat);
A4=kron(Mod.DN(:,nw+1:nw+ny),eye(nw))*kron(eye(ny),uw')*Mod.Hyww(:,1:nx)*Mod.Hxss;
Hy_wzz=zeros(ny*nw*nz,nz);
nyw=ny*nw;
for ii=1:nyw
    Hy_wzz((ii-1)*nz+(1:nz),:)=Hywww((ii-1)*nw+(nx+1:nw),nx+1:nw);
end
A5=kron(Mod.DN(:,nw+1:nw+ny),eye(nw))*tracem(kron(eye(nw*ny),smat)*kron(kron(eye(ny),uw'),eye(nz))*Hy_wzz);
bigmat=kron(Mod.DN(:,1+na+nw:2*na),eye(nw))+kron(Mod.DN(:,nw+1:na),eye(nw))*kron(eye(ny),uw');
bigmat=[bigmat,kron(Mod.DN(:,1:nx),eye(nw))+kron(Mod.DN(:,nw+1:na),eye(nw))*kron(Mod.Hyw(:,1:nx),eye(nw))];
VecEG=linsolve(bigmat,-(A1+A2+A3+A4+A5));
Hyssw=VecEG(1:ny*nw);
Hxssw=VecEG(1+ny*nw:(nx+ny)*nw);

% Hx_sss and Hy_sss 
if isempty(Mod.Skew); Hxsss=zeros(nx,1); Hysss=zeros(ny,1); return; end

Hy_zzt=zeros(ny,nz^2);
Hy_zzz=zeros(ny*nz^2,nz);
for i=1:ny
    Hy_zzt(i,:)=reshape(Mod.Hyww((i-1)*nw+1+nx:i*nw,1+nx:nw),1,[]);
    Hy_zzz((i-1)*nz*nz+1:i*nz*nz,:)=Hy_wzz((i-1)*nw*nz+1+nx*nz:i*nw*nz,:);
end
A1=tracem(kron(kron(eye(nx+ny),nmat'),nmat')*Mod.TN*nmat*Mod.Skew);
A2=3*tracem(kron(eye(nx+ny),[zeros(nx+nz,nz*nz);Hy_zzt;zeros(nw+ny,nz*nz)]')*Mod.HN*nmat*Mod.Skew);
A3=tracem(kron(eye(ny),Mod.Skew)*Hy_zzz);
A3=Mod.DN(:,nw+1:na)*A3;
bigmat=[(Mod.DN(:,1+nw:na)+Mod.DN(:,1+na+nw:2*na)), (Mod.DN(:,1:nx)+Mod.DN(:,1+nw:na)*Mod.Hyw(:,1:nx))];
VecEG=linsolve(bigmat,-(A1+A2+A3));
Hysss=VecEG(1:ny);
Hxsss=VecEG(1+ny:nx+ny);


return;
end


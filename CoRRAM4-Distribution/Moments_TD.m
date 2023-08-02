function Corr = Moments_TD( Mod)
% Computes the correlation matrix between the variables of the canonical DSGE model from the first-order solution

%{
    Copyright: Alfred Mauﬂner

    Purpose: Computes the correlation matrices at lags k=0, 1, ... maxlag between
             the variables s(t):=(x_{1t}, ..., x_{nx,t}, y_{1t}, ..., y_{ny,t}, z_{1t}, ..., z_{nz,t})
             and s(t-k) of the canonical DSGE model. See equation (4.13) of Heer and Mauﬂner,
             Dynamic General Equilibrium Modeling, 3rd. Edition.

    Revision history:

        5 April  2019, first version
       31 May    2019, ordering changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))
       20 August 2019, optionally: solve the Lyapunov equation via the vec operator

    Input:  Mod, an instance of the DSGE class

    Output: Corr, (nx+ny+nz) by (nx+ny+nz) by maxlag+1 array. Corr(i,j,k) stores
            the correlation coefficient between variable i at t and variable j at t-k.

            Corr(i,i,1) stores the standard deviations

 %}

% for simplicity
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nw=nx+nz;
kbar=Mod.maxlag;
nvar=nx+ny+nz;

% initialize
Corr=zeros(nvar,nvar,1+kbar);

% compute Gamma^w_0
Sigma_nu=[zeros(nx,nw);[zeros(nz,nx), Mod.Omega*Mod.Omega']];
Amat=[Mod.Hxw;[zeros(nz,nx), Mod.Rho]];
if Mod.Flags.Sylvester
    Gamma0=Lyapunov(Amat,Sigma_nu);
else
    Mmat=eye(nw^2)-kron(Amat,Amat);
    Gamma0=linsolve(Mmat,Sigma_nu(:));
    Gamma0=reshape(Gamma0,nw,nw);
end

% compute Gamma^s_k
M=[[eye(nx), zeros(nx,nz)];Mod.Hyw;[zeros(nz,nx),eye(nz)]];
if ~Mod.Flags.Log % the model is formulated in levels so we must convert to relative deviations
    %w=[Mod.VarVal(1:nx);Mod.VarVal(1+nx+nz:nvar);Mod.VarVal(1+nx:nx+nz)];
    w=Mod.VarVal;
    w(w==0)=1;
    M=diag(1./w)*M;
end    
Cov=M*Gamma0*M';
Std=sqrt(diag(Cov));
D=diag(1./Std);
Corr(:,:,1)=D*Cov*D;
for ii=1:nvar; Corr(ii,ii,1)=Std(ii); end    
for k=1:kbar
    Cov=M*((Amat^k)*Gamma0)*M';
    Corr(:,:,1+k)=D*Cov*D;
end

return;

end


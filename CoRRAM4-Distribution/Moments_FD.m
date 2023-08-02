function Corr = Moments_FD(Mod)
% Computes the correlation matrix between the variables of the canonical DSGE model from the first-order solution in the
% frequency domain

%{
    Copyright: Alfred Maußner

    Purpose: Computes the correlation matrices at lags k=0, 1, ... maxlag between
             the variables s(t):=(x_{1t}, ..., x_{nx,t}, y_{1t}, ..., y_{ny,t}, z_{1t}, ..., z_{nz,t})
             and s(t-k) of the canonical DSGE model. See Section 4.3.4 of Heer and Maußner,
             Dynamic General Equilibrium Modeling, 3rd. Edition.

    Revision history:

        8 April 2019, first version
       31 May   2019, ordering of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))


    Input:  Mod, an instance of the DSGE class

    Output: Corr, (nx+ny+nz) by (nx+ny+nz) by maxlag+1 array. Corr(i,j,k) stores
            the correlation coefficient between variable i at t and variable j at t-k.

            Corr(i,i,1) stores the standard deviations

    Remarks: 

     If Mod.HPL is different from zero, the function returns moments of HP-filtered  time series.

 %}

% for simplicity
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
kbar=Mod.maxlag;
nvar=nx+ny+nz;

% initialize
nodes=256;
freq=0:(2*pi/nodes):(2*pi*(1-0.5/nodes));
LAMBDA=Mod.hpl;
if Mod.hpl>0;
    hp1= 4*LAMBDA*(1 - cos(freq)).^2 ./ (1 + 4*LAMBDA*(1 - cos(freq)).^2);
    hp=  4*Mod.hpl*(1-cos(freq)).^2 ./ (1+4*Mod.hpl*(1-cos(freq)).^2); % gain of the HP-filter, see eq. (4.??)
else
    hp=ones(1,nodes);
end;
Corr=zeros(nvar,nvar,1+kbar);
Sigma=Mod.Omega*Mod.Omega';
Mmat=[[eye(nx), zeros(nx,nz)];Mod.Hyw;[zeros(nz,nx),eye(nz)]];
if ~Mod.Flags.Log; % the model is formulated in levels so we must convert to relative deviations
    %w=[Mod.VarVal(1:nx);Mod.VarVal(1+nx+nz:nvar);Mod.VarVal(1+nx:nx+nz)];
    w=Mod.VarVal;
    w(w==0)=1;
    Mmat=diag(1./w)*Mmat;
end;    

SW=zeros(nodes,(nx+nz)^2);
im=sqrt(-1);

% power spectrum of process for w(t+1)
for k=1:nodes;
            z = exp(im*freq(k));
           zi = exp(-im*freq(k));
          tmp = [inv(eye(nx)-Mod.Hxw(:,1:nx)*zi)*Mod.Hxw(:,1+nx:nx+nz)*zi;eye(nz)]* ...
                inv(eye(nz)-Mod.Rho*zi)*Sigma*inv(eye(nz)-Mod.Rho'*z)*...
                [Mod.Hxw(:,1+nx:nx+nz)'*z*inv(eye(nx)-Mod.Hxw(:,1:nx)'*z), eye(nz)];
  SW(k,:) = (hp(k)^2)*conj((tmp(:))');
end;
    
MomW=real(ifft(SW));
Cov=Mmat*reshape(MomW(1,:),nx+nz,nx+nz)*Mmat';
Std=sqrt(diag(Cov));
D=diag(1./Std);

for k=1:Mod.maxlag+1;
    tmp=Mmat*reshape(MomW(k,:),nx+nz,nx+nz)*Mmat';
    Corr(:,:,k)=D*tmp*D;
end;

for k=1:nvar;
    Corr(k,k,1)=Std(k);
end;

return;

end


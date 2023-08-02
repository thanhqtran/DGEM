function [FN,DS,DN,HS,HN,TS,TN,rc] = GetFDHT_SD( Mod )
% Compute the Jacobian and if desired the Hessian and the matrix of third-order derivatives of the canonical DSGE model
%{
    Copyright:  Alfred Maußner    

    Purpose:  This function uses the Matlab symbolic toolbox to
              compute numerical derivatives of the system of equations that
              defines an instance of the canonical DSGE model (see SolveModel.m for
              and DSGE.m) from the symbolic expressions for these matrices
              derived the the jacobian command.

    Revision history:          
           
    8 March   2019, first version (from the former script GetDHT)
   15 March   2019, computation of Hessian and third-order matrix added
   24 May     2019, ordering of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))
   21 June    2019, 

    Input: Mod, an instance of the DSGE class
    Output: DS, symbolic Jacobian matrix
            DN, numeric Jacobian matrix
            HS, symbolic Hessian matrix
            HN, numeric Hessian matrix
            TS, symbolic matrix of third-order derivatives
            TN, numeric matrix of third-order derivatives
            rc, return code that defines errors

    
%}

% initialize
HS=[];
HN=[];
TS=[];
TN=[];
rc=0;
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nxy=nx+ny;
na=nxy+nz;

ind=[1:nx,1+nxy:na,1+nx:nxy,1+na:na+nx,1+na+nxy:2*na,1+na+nx:na+nxy];

% define file handle (existence is checked by invoking the DSGE class)
Sys=str2func(Mod.Equations);

% symbolic Jacobian
[FS, x]=Sys();
xn=x(ind); % reordering to [x',z',y',x,z,y]

if ~Mod.Flags.LoadDHT || isempty(Mod.DS)
    DS=jacobian(FS,xn);
    DS=simplify(DS);
else
    DS=Mod.DS;
end

% numeric evaluation
% first: parameters
npar=length(Mod.ParVal);
for is = 1:npar
    eval([Mod.ParSym{is},' = Mod.ParVal(is);']);
end
clear is npar;
% second: variables
nvar=Mod.nx+Mod.nz+Mod.ny;
for is = 1:nvar
    str1=strcat(Mod.VarSym{is},'1');
    str2=strcat(Mod.VarSym{is},'2');    
    eval([str1,'=Mod.VarVal(is);']);
    eval([str2,'=Mod.VarVal(is);']);
end

% third evaluate system of equations at the stationary solution
FN=eval(FS);
if max(abs(FN))>Mod.etol % check stationary solution
    DN=[];
    rc=2;
    return;
end

% fourth: numeric Jacobian
DN=eval(DS);

% Check for invalid Jacobian
if any(any(isnan(DN))) || any(any(isinf(DN)))
    rc=3;
    return;
end

if Mod.order==1; return; end    

% Hessian
if ~Mod.Flags.LoadDHT || isempty(Mod.HS)
    HS=jacobian(DS(:),xn);
    HS=simplify(HS);
else
    HS=Mod.HS;
end
HN=eval(HS);

% arrange
nxy=Mod.nx+Mod.ny;
na=nxy+Mod.nz;
n_m=2*na;  % number of rows of the HN
imat1=zeros(n_m,nxy);
for ii=1:n_m
    imat1(ii,:)=(ii-1)*nxy+1:ii*nxy;
end
HN=HN(imat1,:);

if Mod.order<3; return; end

% Third-order matrix
if ~Mod.Flags.LoadDHT || isempty(Mod.TS)
    TS=jacobian(HS(:),xn);
    TS=simplify(TS);
else
    TS=Mod.TS;
end
TN=eval(TS);

% arrange 
m2=n_m^2;
imat1=zeros(m2,nxy);
for ii=1:m2
    imat1(ii,:)=((ii-1)*nxy+1:ii*nxy);
end
TN=TN(imat1(:),:);

return;

end
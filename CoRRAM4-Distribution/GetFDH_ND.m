function [FN,D,H,rc] = GetFDH_ND( Mod )
% Compute the Jacobian and if desired the Hessian matrix of the canonical DSGE model
%{
    Copyright:  Alfred Maußner    

    Purpose:  This function uses the functions CDJac and CDHesse to
              compute numerical derivatives of the system of equations that
              defines an instance of the canonical DSGE model (see DSGE.m)

    Revision history:          
           
    24 May 2019, first version (from GetFDH_ND)

    Input: Mod, an instance of the DSGE class

    Output: FN, the system of equations evaluated at the stationary solution
             D, the Jacobian matrix
             H, the Hessian matrix
           r c, return code that defines errors

    Revision history:

        2 August    2019, first version
       26 September 2019, bug in line 54 (call of CDJacRE) fixed

%}

% intialize
D=[];
H=[];
rc=0;
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nxy=nx+ny;    
na=nxy+nz;
ind=[1:nx,1+nxy:na,1+nx:nxy,1+na:na+nx,1+na+nxy:2*na,1+na+nx:na+nxy];

% define file handle (existence is checked by invoking the DSGE class)
Sys=str2func(Mod.Equations);    
eq_n=0;    
Sys2=@(x)Sys(x,eq_n,Mod.ParVal);
FN=Sys2([Mod.VarVal;Mod.VarVal]);
if max(abs(FN))>Mod.etol % check stationary solution
    rc=2;
    return;
end
% compute the jacobian matrix of the model
if Mod.Flags.RE
    D=CDJacRE(Sys2,[Mod.VarVal;Mod.VarVal]);
else
    D=CDJac(Sys2,[Mod.VarVal;Mod.VarVal]);
end
D=D(:,ind);

% Check for invalid Jacobian
if any(any(isnan(D))) || any(any(isinf(D)))
    rc=3;
    return;
end
if Mod.order==1
    return;
else
    H=zeros(nxy*2*na,2*na);
    for eqn=1:nxy
        Sys2=@(x)Sys(x,eqn,Mod.ParVal);
        if Mod.Flags.RE
            temp=CDHesseRE(Sys2,[Mod.VarVal;Mod.VarVal]);
        else
            temp=CDHesse(Sys2,[Mod.VarVal;Mod.VarVal]);
        end
        H(1+(eqn-1)*2*na:eqn*2*na,:)=temp(ind,ind);    
    end
end    
if Mod.order<3; return; end
end


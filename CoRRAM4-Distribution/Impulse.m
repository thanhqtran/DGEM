function IRF = Impulse(Mod)
%Computes impulse respones for an instance of the canonical DSGE model in Mod

%{
    Copyright: Alfred Mauﬂner
    
    Purpose:  this function computes impulse respones for all variables
              of an instance of the canonical DSGE model specified in Mod.

              The impulse responses are computed only from the linear part of the solution
              and are perentage deviations of each variable relative to the unshocked 
              path.

              See Chapter 4 of Heer and Mauﬂner, Dynamic General Equilibrium Modeling, 3rd
                  edition. 

    Revision history:

              3 April 2018 (from previous versions in SimulateModel.m)
             31 May   2019, order of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))

    Input:    Mod, an instance of the DSGE class

    Output:   IRF, nobs by nx+ny+1 by bz array. Element IRF(t,i,s) stores
              the percentage deviation of variable i at time t in response
              to a one-time shock at time t=2 to variable s. In position i=nx+ny+1 the
              array stores the time path of the shocked variable s=1, 2, ..., nz,
              where nz is the number of shocks, nx the number of endogenous state
              variables, and ny the number of jump variables.

%}

% initialize
nobs=Mod.inobs; % number of periods
IRF=zeros(nobs,Mod.nx+Mod.ny+1,Mod.nz); % array with solutions
eunit=eye(Mod.nz);
if ~Mod.Flags.Log % store stationary values and multipliers
    v0=Mod.VarVal(1:Mod.nx+Mod.ny);
    v0(v0==0)=1; % replace values equal to zero with ones
    m0=ones(1,Mod.nx+Mod.ny)*100;
    m0(v0==0)=1;
    m0=diag(m0); % values with zero stationary solution remain in absolute deviations
else
    m0=diag(ones(1,Mod.nx+Mod.ny)*100); % log model, zero stationary value not allowed
end
if Mod.Flags.Ds % xi
    m1=diag(Mod.VarXi(1:Mod.nx+Mod.ny));
end

for ns=1:Mod.nz % outer loop over shocks
    % initialize
    zt=zeros(Mod.nz,1+nobs);
    xt=zeros(Mod.nx,1+nobs);
    yt=zeros(Mod.ny,nobs);
    zt(:,2)=Mod.Omega*eunit(:,ns);
    for tt=2:nobs
        zt(:,tt+1)=Mod.Rho*zt(:,tt);
        xt(:,tt+1)=Mod.Hxw*[xt(:,tt);zt(:,tt)];
        yt(:,tt)  =Mod.Hyw*[xt(:,tt);zt(:,tt)];
    end
    if ~Mod.Flags.Log % compute percentage deviations
        xyt=bsxfun(@rdivide,[xt(:,1:nobs);yt],v0);
    else
        xyt=[xt(:,1:nobs);yt];
    end
    if Mod.Flags.Ds % add cumulative sum of growth factor
        ahut=repmat(cumsum([0 xyt(1+Mod.nx,1:nobs-1)]),Mod.nx+Mod.ny,1);
        xyt=xyt+m1*ahut;
    end
    IRF(:,:,ns)=[(m0*xyt)',100*zt(ns,1:nobs)'];
        
end


end


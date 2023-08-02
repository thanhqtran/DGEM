function smat = Simulate_NP(Mod,EpsMat,s,PF)
%Simulates a time path for a DSGE model without pruning
%{
    Alfred Mauﬂner

    Purpose: simulate a time path for a DSGE model.

    Input:   Mod, an instance of the DSGE class, see DSGE.m

             EpsMat, nz by nobs+1 matrix with random numbers,
                     where nz is the number of shocks

             s, scalar, number of simuation

    Output:  smat, nobs by nvar matrix, each column stores the simulated time series
             for one of the nvar variables of the model. The order of the columns is
             x(1) - x(nx), y(1), ..., y(ny), z(1), ..., z(nz), where x are the model's
             endogenous states, y the model's jump variables, and z the shocks.

    Revision history:

        12 April    2019, first version (from a previous version of SimulateModel.m)
        18 April    2019, nobs periods for burn in and nobs for the estimation of second moments
        31 May      2019, ordering of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))
        06 August   2019, difference stationary growth is no longer treaded here but in Filtering.m,
                          consistent with Simulate_P2.m and Simulate_P3.m
        09 December 2020, policy function for simulation added as additional argument to allow for
                          non-perturbations solutions
%}

% for simplicity
nobs=Mod.nobs+Mod.burnin;
smat=[];

% initialize
xt=zeros(Mod.nx,nobs+1);
yt=zeros(Mod.ny,nobs);
zt=zeros(Mod.nz,nobs+1);

if Mod.Flags.CheckBounds
    x_l=Mod.VarBnd(1:Mod.nx,1);
    x_u=Mod.VarBnd(1:Mod.nx,2);
    y_l=Mod.VarBnd(1+Mod.nx:Mod.nx+Mod.ny,1);
    y_u=Mod.VarBnd(1+Mod.nx:Mod.nx+Mod.ny,2);
    fstring='%3i %10.5f %10.5f %10.5f\n';
end

zt(:,1)=Mod.Omega*EpsMat(:,1);
xt(:,1)=Mod.VarVal(1:Mod.nx);

if ~Mod.Flags.CheckBounds
    % iterate without checking of bounds
    for t=1:nobs
        zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*EpsMat(:,t+1);
        [xt(:,t+1),yt(:,t)]=PF([xt(:,t);zt(:,t)]);
        if any(isnan(xt(:,t+1))); error('Nan in xt at simulation s=%u and period t=%u',s,t); end
    end
else
    for t=1:nobs
        zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*EpsMat(:,t+1);
        %[xt(:,t+1),yt(:,t)]=Mod.PFXY([xt(:,t);zt(:,t)]);
        [xt(:,t+1),yt(:,t)]=PF([xt(:,t);zt(:,t)]);
        if ~WithinBounds(xt(:,t+1),yt(:,t),t,s); return; end
    end
end

% compute deviations (if the model is in logs these are equal to relative
% deviations
xt=xt(:,1+Mod.burnin:Mod.nobs+Mod.burnin)-Mod.VarVal(1:Mod.nx);
yt=yt(:,1+Mod.burnin:Mod.nobs+Mod.burnin)-Mod.VarVal(1+Mod.nx:Mod.nx+Mod.ny);

if ~Mod.Flags.Log % compute relative deviations, if the model in in levels
    v0=Mod.VarVal(1:Mod.nx+Mod.ny);
    v0(v0==0)=1; % replace values equal to zero with ones
    smat=bsxfun(@rdivide,[xt;yt],v0);  
else
    smat=[xt;yt];
end
smat=[smat;zt(:,1+Mod.burnin:Mod.nobs+Mod.burnin)]';

return;    

    function ok=WithinBounds(x1,y1,t,s)
        ind1=find(x1<x_l|x1>x_u);
        if ~isempty(ind1) 
            tmp=[x_l,x1,x_u];
            fprintf(Mod.lfh,'%s %3i %s %3i %s\n','Simulation s=',s,' time t=',t,' state variable(s) out of bounds');            
            fprintf(Mod.lfh,fstring,[ind1,tmp(ind1,:)]');
            fprintf(Mod.lfh,'%s\n',' ');
            ok=false;
            return;
        else
            ind1=find(y1<y_l|y1>y_u);
            if ~isempty(ind1)
                tmp=[y_l,y1,y_u];
                fprintf(Mod.lfh,'%s %3i %s %3i %s\n','Simulation s=',s,' time t=',t,' control variable(s) out of bounds');
                fprintf(Mod.lfh,fstring,[ind1,tmp(ind1,:)]');
                fprintf(Mod.lfh,'%s\n',' ');
                ok=false;
                return;
            end
        end
        ok=true;
        return;
    end
end


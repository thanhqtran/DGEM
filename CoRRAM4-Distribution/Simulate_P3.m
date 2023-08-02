function smat = Simulate_P3(Mod,EpsMat,s)
%Simulates a time path for a DSGE model from the thrid-order solution with pruning
%{
    Alfred Mauﬂner

    Purpose: simulate a time path for a DSGE model from the third-order, pruned solution.

    Input:   Mod, an instance of the DSGE class, see DSGE.m

             EpsMat, nz by nobs+1 matrix with random numbers,
                     where nz is the number of shocks

             s, scalar, number of simuation

    Output:  smat, nobs by nvar matrix, each column stores the simulated time series
             for one of the nvar variables of the model. The order of the columns is
             x(1) - x(nx), y(1), ..., y(ny), z(1), ..., z(nz), where x are the model's
             endogenous states, y the model's jump variables, and z the shocks.

    Revision history:

        25 April 2019, first version (from Simulate_P2.m)
        31 May   2019, ordering of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))
        06 August 2019, shock added to smat
%}

% for simplicity
nobs=Mod.nobs+Mod.burnin;
smat=[];

% initialize
xft=zeros(Mod.nx,nobs+1);
xst=xft;
xrt=xft;

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

if ~Mod.Flags.CheckBounds
    % iterate without checking of bounds
    for t=1:nobs
        zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*EpsMat(:,t+1);
        xft(:,t+1)=Mod.Hxw*[xft(:,t);zt(:,t)];
        xst(:,t+1)=Mod.Hxw*[xst(:,t);zt(:,t)] + 0.5*kron(eye(Mod.nx),[xft(:,t);zt(:,t)]')*Mod.Hxww*[xft(:,t);zt(:,t)]+0.5*Mod.Hxss;
        xrt(:,t+1)=Mod.Hxw*[xrt(:,t);zt(:,t)] + kron(eye(Mod.nx),[xft(:,t);zt(:,t)]')*Mod.Hxww*xst(:,t)*[xft(:,t);zt(:,t)] ...
                   + (1/6)*Mod.Hxsss + 0.5*Hxssw*[xft(:,t);zt(:,t)] ...
                   + (1/6)*kron(kron(eye(Mod.nx),[xft(:,t);zt(:,t)]'),[xft(:,t);zt(:,t)]')*Mod.Hxwww*[xft(:,t);zt(:,t)];
        yt(:,t)=Mod.Hyw*[(xft(:,t)+xst(:,t))+xrt(:,t);zt(:,t)] ...
                + 0.5*kron(eye(Mod.ny),[xft(:,t);zt(:,t)]')*Mod.Hyww*([xft(:,t);zt(:,t)]+2*[xst(:,t);zt(:,t)]) ...
                + 0.5*Mod.Hyss + (1/6)*Mod.Hysss+0.5*Mod.Hyssw*[xft(:,t);zt(:,t)] ...
                + (1/6)*kron(kron(eye(Mod.ny),[xft(:,t);zt(:,t)]'),[xft(:,t);zt(:,t)]')*Mod.Hywww*[xft(:,t);zt(:,t)];
    end
else
    xt=xft;
    xstar=Mod.VarVal(1:Mod.nx);
    ystar=Mod.VarVal(1+Mod.nx:Mod.nx+Mod.ny);
    for t=1:nobs
        zt(:,t+1)=Mod.Rho*zt(:,t)+Mod.Omega*EpsMat(:,t+1);
        xft(:,t+1)=Mod.Hxw*[xft(:,t);zt(:,t)];
        xst(:,t+1)=Mod.Hxw*[xst(:,t);zt(:,t)] + 0.5*kron(eye(Mod.nx),[xft(:,t);zt(:,t)]')*Mod.Hxww*[xft(:,t);zt(:,t)]+0.5*Mod.Hxss;
        xrt(:,t+1)=Mod.Hxw*[xrt(:,t);zt(:,t)] + kron(eye(Mod.nx),[xft(:,t);zt(:,t)]')*Mod.Hxww*xst(:,t)*[xft(:,t);zt(:,t)] ...
                   + (1/6)*Mod.Hxsss + 0.5*Hxssw*[xft(:,t);zt(:,t)] ...
                   + (1/6)*kron(kron(eye(Mod.nx),[xft(:,t);zt(:,t)]'),[xft(:,t);zt(:,t)]')*Mod.Hxwww*[xft(:,t);zt(:,t)];
        xt(:,t+1) = xft(:,t+1)+xst(:,t+1)+xrt(:,t+1);
               yt(:,t)=Mod.Hyw*[xt(:,t);zt(:,t)] ...
                + 0.5*kron(eye(Mod.ny),[xft(:,t);zt(:,t)]')*Mod.Hyww*([xft(:,t);zt(:,t)]+2*[xst(:,t);zt(:,t)]) ...
                + 0.5*Mod.Hyss + (1/6)*Mod.Hysss+0.5*Mod.Hyssw*[xft(:,t);zt(:,t)] ...
                + (1/6)*kron(kron(eye(Mod.ny),[xft(:,t);zt(:,t)]'),[xft(:,t);zt(:,t)]')*Mod.Hywww*[xft(:,t);zt(:,t)];
        if ~WithinBounds(xstar+xt(:,t+1),ystar+yt(:,t),t,s); return; end
    end
end

xrt=xrt(:,1+Mod.burnin:Mod.nobs+Mod.burnin)-Mod.VarVal(1:Mod.nx);
yt=yt(:,1+Mod.burnin:Mod.nobs+Mod.burnin)-Mod.VarVal(1+Mod.nx:Mod.nx+Mod.ny);

% compute relative deviations, if the model is in levels
if ~Mod.Flags.Log 
    v0=Mod.VarVal(1:Mod.nx+Mod.ny);
    v0(v0==0)=1; % replace values equal to zero with ones
    %smat=bsxfun(@rdivide,[xst;yt],v0);  
    smat=[xrt;yt]./v0;
else
    smat=[xrt;yt];
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


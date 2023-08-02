function xf = Filtering(Mod,x)
%Returns the cyclical component of the columns of x computed from various filters
%{
    Alfred Mauﬂner

    Purpose: Computes for each column of the matrix x its cyclical component obtained
             from passing it through different filters.
 
             The data in x are assumed to be relative deviations of the model's variables
             from their value at the model's deterministic steady-state.

             In a first step we compute the time paths for all variables
             depending on the setting of the flag Mod.Flags.Ds. In the second step 
             we path these through various filters.

     Input:  Mod, an instance of the DSGE class (see, DSGE.m)
             x, nobs by nvar matrix. The first Mod.nx columns pertain
                to simulated time series for the model's endogenous state variables,
                the next Mod.ny pertain to simulated time series for the model's jump
                variables and the remaining Mod.nz columns to the model's shocks.

     Output: xf, nobs by nvar matrix, the filtered series

     Revisions: 06 August 2019, bugs in data preparation for difference stationary models corrected
             
%}
% get dimensions
[nobs,nvar]=size(x);
nxy=nvar-Mod.nz;

if ~Mod.Flags.Ds
    switch Mod.Filter
        case 'HP'
            xf=HPF(x,Mod.hpl);
        case 'Linear'
            xf=x;
        case 'H'
            xf=HamiltonFilter(x);
            %m1=diag([Mod.VarXi(1:Mod.nx);Mod.VarXi(1+Mod.nx+Mod.nz:Mod.nx+Mod.nz+Mod.ny)]);
            %tmp=m1*repmat((1:1:nobs)*log(Mod.ParVal(1)),nxy,1)+x(:,1:nxy)'; % growth factor is first element in the parameter vector    
            %xf=HamiltonFilter([tmp',x(:,1+nxy:nvar)]);
    end
else
    m1=diag(Mod.VarXi(1:Mod.nx+Mod.ny));
    tmp=x(:,1:nxy)';
    ahut=repmat(cumsum([0 tmp(1+Mod.nx,1:nobs-1)]),Mod.nx+Mod.ny,1);
    tmp=tmp+m1*ahut;
    tmp=[tmp',x(:,1+nxy:nvar)];
    switch Mod.Filter
        case 'HP'
            xf=HPF(tmp,Mod.hpl);
        case 'Linear'
            xf=tmp;
        case 'H'
            at=log(Mod.ParVal(1))+x(:,1+Mod.nx); % deviations of the growth factor are first element in the vector of jump variables
            at=repmat(cumsum([0 at(1:nobs-1)]),nxy,1);
            tmp=m1*at+x(:,1:nxy)';
            tmp=[tmp',x(:,1+nxy:nvar)];
            xf=HamiltonFilter(tmp);
    end
end


return

end
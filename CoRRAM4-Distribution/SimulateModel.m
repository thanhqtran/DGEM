function Mod = SimulateModel(Mod,varargin)
% Computes array with impulse reponses and arrays with second moments for a
% DSGE model

%{
    Alfred Maußner

    Purose: computes an array with impulse respones for each shock
            and arrays with second moments for an instance of the DSGE class (see DSGE.m)

    Input: Mod, an instance of the DSGE class

    Output: Mod, the same instance with solutions added to the repective class properties.

    Revison history:

        12 April    2019, first version
        18 April    2019, burnin phase added
         9 May      2019, IRFs and second moments are optional
        09 December 2020, policy function from another solution procedure as optional argument added 

%}

% impulse response
if Mod.Flags.IRFs
    Mod.IRF=Impulse(Mod);
end

% second moments
if ~Mod.Flags.Moments; return; end
if Mod.Flags.Simulate
    % check optional input
    if nargin==1
        PF=@(v)PFXY(Mod,v); % employ policy function from the perturbation solution
    else
        PF=varargin{1}; % employ policy function provided as optional argument
    end
    rng(Mod.seed);
    epsMat=randn(Mod.nz,Mod.nobs+1+Mod.burnin,Mod.nofs);
    nvar=Mod.nx+Mod.nz+Mod.ny;
    Mod.Corr_m=zeros(nvar,nvar,Mod.maxlag+1);
    if Mod.Flags.Confidence
        CorMat=zeros(nvar,nvar,1+Mod.maxlag,Mod.nofs);
    end
    Mod.valsim=0;
    for s=1:Mod.nofs
        switch Mod.Prune
            case 0
                tmp=Simulate_NP(Mod,epsMat(:,:,s),s,PF); % 9 Dec 2020: PF as additional input added  
            case 2
                if isempty(Mod.Hxww); error('The model has not been solved to the second-order'); end
                tmp=Simulate_P2(Mod,epsMat(:,:,s),s);
            case 3
                if isempty(Mod.Hxwww); error('The model has not been solved to the thrid-order'); end
                tmp=Simulate_P3(Mod,epsMat(:,:,s),s);
        end
        if ~isempty(tmp)
            if ~strcmp('None',Mod.Filter)
                tmp=Filtering(Mod,tmp);
            end 
            Mod.valsim=Mod.valsim+1;
            tmp=CCR(tmp,Mod.maxlag);
            Mod.Corr_m=Mod.Corr_m+tmp;
            if Mod.Flags.Confidence
                CorMat(:,:,:,Mod.valsim)=tmp;
            end
        end
    end
    if Mod.valsim>0
        Mod.Corr_m=Mod.Corr_m/Mod.valsim;
        if Mod.Flags.Confidence
            tmp=sort(CorMat(:,:,:,1:Mod.valsim),4);
            imin=floor(Mod.nofs*0.025);
            imax=round(Mod.nofs*0.975);
            Mod.Corr_l=reshape(tmp(:,:,:,imin),nvar,nvar,1+Mod.maxlag,[]);
            Mod.Corr_u=reshape(tmp(:,:,:,imax),nvar,nvar,1+Mod.maxlag,[]);
        end
    else
        Mod.Corr_l=[]; Mod.Corr_u=[]; Mod.Corr_m=[];
    end
else
    if strcmp('HP',Mod.Filter)
        Mod.Corr_m=Moments_FD(Mod);
    else
        Mod.Corr_m=Moments_TD(Mod);
    end
end
return;

end


classdef DSGE
    % Definition of the canonical dynamic stochastic general equilibrium model
    %{
        Copyright: Alfred Maußner
    
        Purpose:  this class summarizes the properties of the canonical DSGE model
                  of Chapter 2 of Heer and Maußner, Dynamic General Equilibrium Modeling, 3rd
                  edition. In addition, it also holds the mode's solution and simulation results.
    
        Revisions: 10 January   2017, first version
                   11 May       2017, non zero means added
                   19 February  2019, wstar added to hold stationary solution
                   12 March     2019, fundamentally revised according to CoRRAM version 4 for Gauss
                   27 March     2019, 3rd-order added, Flags.numeric deleted, information about derivatives is stored in 
                                      Derivatives
                   25 April     2019, Filter names added as property
                   31 May       2019, Flag for Richardson extrapolation added
                   02 August    2019, stationary solution in w now an additional optional argument
                   13 August    2019, flag for write to Excel spreadsheet added
                   25 September 2019, vector valued policy functions up to order=2 added
                   30 September 2019, Loadpath added to path, if not already in file list
                   21 November  2019, in PFXY_v check for nx,nz>0
                   28 November  2020, in PFXY_v take care of nx=0 and nz=1 as well as of nx=1 and nz=0
                   12 February  2021, in PFXY, transpose if input is a row instead of a column vector
                   16 February  2022, Flag.Append added and initilized at false
    
    %}
    
    properties
        nx;             % number of endogenous state variables
        ny;             % number of not predetermined variables
        nz;             % number of pure shocks
        nu;             % number of static equations
        Rho=[];         % autocorrelation matrix of shocks        
        Omega=[];       % standard deviations of innovations        
        Skew=[];        % nz by nz^2 skewness matrix
        Mu=[];          % nz by one vector, the non-zero means of the innovations
        Muss=[];        % nz by one vector of second order derivatives of mean vector of innovations                
        Equations='';   % string the name of the file that defines the model's equations
        Var=struct('Name','','Symbol','','Type','','Pos',0,'Print',0,'Corr',0,'Rel',0,'Star',0,'Xi',0,'Plotno',0,'Bound',[0 0]);    
        ParVal;         % place to pass parameters to the equations of the model
        ParSym;         % place to pass parameter names to the equations of the model
        VarVal=[];      % stores stationary solution        
        VarBnd=[];      % stores upper and lower bounds for simulation
        VarSym;         % cell array for variable symbols
        VarXi=[];       % stores the scaling factors
        valsim=[];      % stores the number of valid simulations
        FN=[];          % the solution of the system of equations at the stationary solution
        DN=[];          % numeric Jacobian matrix
        DS=[];          % symbolic Jacobian matrix
        HN=[];          % numeric Hessian
        HS=[];          % symbolic Hessian
        TS=[];          % symbolic matrix of third-order derivatives
        TN=[];          % numeric matrix of third-order derivatives
        Hxw=[];         % matrix of first-order coefficients for the endogenous states
        Hyw=[];         % matrix of first-order coefficients for the controls
        Hxww=[];        % matrix of second-order coefficients for the endogenous states
        Hyww=[];        % matrix of second-oder coefficients for controls
        Hxss=[];        % vector of second-order coefficients for the endogenous states wrt the perturbation parameter
        Hyss=[];        % vector of second-order coefficients for controls wrt to the perturbation parameter
        Hxwww=[];       % matrix of third-order coefficients for the endogenous states wrt to all states
        Hywww=[];       % matrix of third-order coefficients for the controls wrt to all states
        Hxssw=[];       % matrix of third-order coefficients for the endogenous states wrt to state dependend risk
        Hyssw=[];       % matrix of third-order coefficients for the controls wrt to state dependend risk
        Hxsss=[];       % vector of third-order coefficients for the endogenous states wrt the perturbation parameter
        Hysss=[];       % vector of third-order coefficients for the control wrt the perturbation parameter
        Corr_m=[];      % array with second moments
        Corr_l=[];      % array with lower bounds
        Corr_u=[];      % array with upper bounds
        IRF=[];         % array with impulse respones      
        Flags=struct('Balance',0,'CheckBounds',0,'Confidence',0,'Ds',0,'Grid',1,'IRFs',1,'LegendBox',1,'Log',0,'LoadDHT',0,'Moments',1,'Reduce',1,'RE',0,'Simulate',1,'Sylvester',1,'Trendline',0,'Excel',0,'Append',0);
        order=1;        % order of solution, 1,2, or 3
        Prune=0;        % pruned system: 0=no pruning, 2=second order, 3=third order
        Filter='HP';    % allowed strings: 'Linear', 'None','HP','H'
        Derivatives='ND'; % string can be: 'AD' for automatic differentiation, 'ND' for numeric and 'SD' for symbolic
        etol=1e-9;      % tolerance used to check whether the model's equations hold at the stationary point
        itol=1e-10;     % tolerance used to check for complex solutions
        inobs=10;       % number of periods for impulse responses
        maxlag=1;       % maximum number of lags considered in the computation of covariances
        burnin=50;      % additional simulations for burn in
        nobs=120;       % number of periods in simulations of the model
        nofs=500;       % number of simulations
        hpl=1600;       % Hodrick-Prescott filter weight
        seed=181191;    % seed for random number generation
        Outfile='Model';% base name of the file that receives tables
        Loadpath=cd;    % path from which Matlab looks for the model
        lfh;            % handle of the logfile
        rc=0;           % error code number
        Messages=cellstr(char('Model does not exist. Check filename and path.',...
                              'Invalid stationary solution', ...
                              'Jacobian includes NaNs and/or Infs', ...
                              'Schur failed', ...
                              'Instable model', ...
                              'Indetermined model'));

    end
    
    methods
        function obj = set.nx(obj,nx)
            if ~isnumeric(nx); error('nx must be non-negative integer value'); else obj.nx=nx; end %#ok<SEPEX>
        end
        function obj = set.ny(obj,ny)
            if ~isnumeric(ny); error('ny must be non-negative integer value'); else obj.ny=ny; end %#ok<SEPEX>
        end
        function obj = set.nz(obj,nz)
            if ~isnumeric(nz); error('nz must be non-negative integer value'); else obj.nz=nz; end %#ok<SEPEX>
        end
        function obj = set.nu(obj,nu)
            if ~isnumeric(nu); error('nu must be non-negative integer value'); else obj.nu=nu; end %#ok<SEPEX>
        end
        function obj = set.Derivatives(obj,str)
            if strcmp(str,'AD') || strcmp(str,'ND') || strcmp(str,'SD')
                obj.Derivatives=str;
            else
                error(strcat('The property Derivatives must be one out of the strings AD, ND and SD. But you set it to ',str));
            end
        end
        
        function rc=ExistsModel(obj)
            str=strcat(obj.Loadpath,'\',obj.Equations,'.m');
            rc=exist(str,'file'); % if true, file exists on Loadpath
            if ~strcmpi(cd,obj.Loadpath); addpath(obj.Loadpath); end % add to path
        end
        
        % this function accepts only a vector v=[x,z] and employs x-xstar to compute the policy functions for
        % x and y. It assumes that min(nx+nz)>=1
        function [x,y]=PFXY(obj,v)
            if min(obj.nx+obj.nz)<1; error('nx and nz cannot be both equal to zero'); end
            if size(v,2)>1; v=v';end
            if obj.nx==0
                vbar=v;
            else
                vbar=[(v(1:obj.nx)-obj.VarVal(1:obj.nx));v(1+obj.nx:obj.nx+obj.nz)];
                x=obj.VarVal(1:obj.nx)+obj.Hxw*vbar;
            end
            if obj.nz==0
                vbar=(v(1:obj.nx)-obj.VarVal(1:obj.nx));
            end
            y=obj.VarVal(1+obj.nx:obj.nx+obj.ny)+obj.Hyw*vbar;
            if obj.order>1.5
                if obj.nx>0
                    x=x+0.5*obj.Hxss+0.5*kron(eye(obj.nx),vbar')*obj.Hxww*vbar;
                end
                y=y+0.5*obj.Hyss+0.5*kron(eye(obj.ny),vbar')*obj.Hyww*vbar;
            end
            if obj.order>2.5                
                if obj.nx>0
                    x=x+(1/6)*obj.Hxsss+0.5*kron(eye(obj.nx),vbar')*obj.Hxssw+(1/6)*kron(kron(eye(obj.nx),vbar'),vbar')*obj.Hxwww*vbar;
                end
                y=y+(1/6)*obj.Hysss+0.5*kron(eye(obj.ny),vbar')*obj.Hyssw+(1/6)*kron(kron(eye(obj.ny),vbar'),vbar')*obj.Hywww*vbar;
            end
            return;
        end
        % this function accepts a (nx+nz) by m matrix vmat, where the first
        % nx rows must contain x, and returns  nx by m and ny by m
        % matrixes with corresponding solutions for the states x and the
        % jump variables y. 
        function[x,y]=PFXY_v(obj,v)
            nx=obj.nx;  %#ok<*PROPLC>
            ny=obj.ny; 
            nz=obj.nz; 
            nw=nx+nz;
            m =size(v,2); 
            % check for proper input
            if min(nx+nz)<1; error('nx and nz cannot be both equal to zero.'); end
            if (nx>0)
                vmat=v-repmat([obj.VarVal(1:nx);zeros(nz,1)],1,m);
                x=repmat(obj.VarVal(1:nx),1,m)+obj.Hxw*vmat;
            else
                vmat=v;
            end
            y=repmat(obj.VarVal(1+nx:nx+ny),1,m)+obj.Hyw*vmat;
            if obj.order>1.5 % for the formula to work, nx+nz must be larger than one so I must take care of this case
                if (nx==1) && (nz==0)
                    x=x+0.5*obj.Hxww*(vmat.^2);
                    y=y+0.5*obj.Hyww*(vmat.^2);
                elseif (nx==0) && (nz==1)
                    y=y+0.5*repmat(obj.Hyss,1,m)+0.5*obj.Hyww*(vmat.^2);
                else
                    x=x+0.5*repmat(obj.Hxss,1,m)+ ...
                    0.5*reshape(dot(reshape(repmat(vmat,nx,1),nw,nx*m),reshape(obj.Hxww*vmat,nw,nx*m)),nx,m);
                    y=y+0.5*repmat(obj.Hyss,1,m)+ ...
                    0.5*reshape(dot(reshape(repmat(vmat,ny,1),nw,ny*m),reshape(obj.Hyww*vmat,nw,ny*m)),ny,m);
                end
            end
            if obj.order>2.5
                error('Not implemented for order 3!');
            end
            return;
        end
        
        function M=DSGE(nx,ny,nz,nu,varargin)
            p=inputParser;
            p.addRequired('nx',@isnumeric);
            p.addRequired('ny',@isnumeric);
            p.addRequired('nz',@isnumeric);
            p.addRequired('nu',@isnumeric);
            p.addOptional('w',@isnumeric);
            p.addOptional('Symbols','',@iscell);
            p.addParameter('Names','',@iscell);
            p.parse(nx,ny,nz,nu,varargin{:});
            nin=length(varargin);
            switch nin
                case 1
                    w=p.Results.w;
                    Symbols=[];
                case 2
                    w=p.Results.w;
                    Symbols=p.Results.Symbols;
                otherwise
                    w=[];
                    Symbols=[];
            end
            Names=p.Results.Names;
                if nargin<4
                    error('Arguments nx, ny, nz, and nu are required');
                elseif nargin>=4
                    M.nx=nx;
                    M.ny=ny;
                    M.nz=nz;
                    M.nu=nu;
                    if ~isempty(w)
                        if length(w)<nx+ny
                            error(strcat('The vector w must have ',num2str(nx+ny,'%4i'),' elements'));
                        end
                        M.VarVal=[w;zeros(nz,1)];
                    end
                    if ~isempty(Symbols); M.VarSym=Symbols; end
                    if (nu>0 && nu<ny); M.Flags.Reduce=true; else M.Flags.Reduce=false; end %#ok<SEPEX>
                    if nx>0
                        for i=1:nx
                            if ~isempty(Names)
                                M.Var(i).Name=char(Names(i));
                            else
                                M.Var(i).Name=strcat('X',num2str(i));
                            end
                            if ~isempty(Symbols)
                                M.Var(i).Symbol=char(Symbols(i));
                            else
                                M.Var(i).Symbol=strcat('X',num2str(i));
                            end
                            M.Var(i).Type='x';
                            M.Var(i).Pos=i;
                            if ~isempty(w); M.Var(i).Star=w(i); end
                        end
                    end
                    if ny>0
                        for i=1:ny
                            if ~isempty(Names)
                                M.Var(i+nx).Name=char(Names(i+nx));
                            else
                                M.Var(i+nx).Name=strcat('Y',num2str(i));
                            end
                            if ~isempty(Symbols)
                                M.Var(i+nx).Symbol=char(Symbols(i+nx));
                            else
                                M.Var(i+nx).Symbol=strcat('Y',num2str(i));
                            end
                            M.Var(i+nx).Type='y';
                            M.Var(i+nx).Pos=i;
                            if ~isempty(w); M.Var(i+nx).Star=w(i+nx); end
                        end
                    end                
                    if nz>0
                        M.Rho=zeros(nz,nz);
                        M.Omega=zeros(nz,nz);
                        M.Skew=zeros(nz,nz*nz);
                        for i=1:nz
                            if ~isempty(Names)
                                M.Var(i+nx+ny).Name=char(Names(i+nx+ny));
                            else
                                M.Var(i+nx+ny).Name=strcat('Z',num2str(i));
                            end
                            if ~isempty(Symbols)
                                M.Var(i+nx+ny).Symbol=char(Symbols(i+nx+ny));
                            else
                                M.Var(i+nx+ny).Symbol=strcat('Z',num2str(i));
                            end
                            M.Var(i+nx+ny).Type='z';
                            M.Var(i+nx+ny).Pos=i;
                            M.Var(i+nx+ny).Star=0;
                            M.Var(i+nx+ny).Plotno=1;
                        end
                    end
                    na=nx+nz+ny;
                    for i=1:na
                        M.Var(i).Plotno=0;
                        M.Var(i).Print=0;
                        M.Var(i).Corr=0;
                        M.Var(i).Rel=0;
                        M.Var(i).Xi=0;
                        M.Var(i).Bound(1)=0; M.Var(i).Bound(2)=0;
                    end
                end
        end

    end
    
end



function  MakeGraphs(Mod,FontName,varargin)
% Plot impulse responses

%{
    Author: Alfred Maußner
    
    Purpose : Plot the impulse responses of an instance of the DSGE class.

    Input   : Mod, and instance of the DSGE class
              FontName, string, the name of the font used to print numbers, labels and legends

    Revision history:
        
        10 May    2019, frist version (from MyPlot), 
        14 August 2019, optional arguments added
    
             
%}

% check whether the model has indeed been simulated
if isempty(Mod.IRF)
    error('There are no IRFs to be plotted');
end

% for simplicity
nx=Mod.nx;
ny=Mod.ny;
nz=Mod.nz;
nvar=nx+ny+nz;

for shock=1:nz
    switch nargin
        case 2 % check the Var structure
            [allNames,ipl]=GetMaxPlotNumber(Mod.Var,shock);
        case {3,4}
            Names=varargin{1};
            ip1=varargin{2};
            nip1=size(ip1,1)-nz;
            ipl=[ip1(1:nip1,:);ip1(nip1+shock,:)];            
            allNames=cell(nvar,1);
            allNames(ipl(1:nip1,1))=Names(1:nip1);
            allNames(nx+ny+1:1:nvar)=Names(nip1+1:nip1+nz);
    end

    % gather information for each figure and plot
    np=max(ipl(:,2));
    panels=1:np;
    
    Pvar=struct('n','','x',[]);
    for jj=1:np
        k1=find(ipl(:,2)==panels(jj));
        if ~isempty(k1)
            k2=ipl(k1,1);
            nk2=length(k2);
            Pvar(jj).n=allNames(k2);            
            if (k2(nk2)==(nx+ny+shock)) % the shock is in the current plot list so plot
                k2(nk2)=nx+ny+1;
            end
            Pvar(jj).x=Mod.IRF(:,k2,shock);
        else
            warning('Plotnumbers should run from 1 to %i but plotnumber %i is not specified!',np,panels(jj));
        end
    end
    figname=allNames{nx+ny+shock};
    figure('Name',figname);

    switch np
            case 1
                Plot1(Pvar,1,1,FontName,10,1.2,'Box',Mod.Flags.LegendBox,'Grid',Mod.Flags.Grid);
            case 2
                Plot1(Pvar,2,1,FontName,10,1.2);
            case {3, 4}
                Plot1(Pvar,2,2,FontName,10,1.2);
            case {5 6}
                Plot1(Pvar,3,2,FontName,9,1.2);
            case {7 8}
                Plot1(Pvar,4,2,FontName,8,1.2);
    end
end

return;
end

% Gather information about the number of panels and the number of lines per panel
function [Names,ipl] = GetMaxPlotNumber(Var,si)

    nvar=length(Var);
    Names=cell(nvar,1);
    ipl=zeros(nvar,2);
    
    for i2=1:nvar
        if Var(i2).Plotno>0
            if Var(i2).Plotno>8
                message=strcat('CoRRAM supports only 8 panels but Var(',num2str(i2),').Plotno=',num2str(Var(i2).Plotno));
                error(message);    
            else
                if (strcmp(Var(i2).Type,'z') && (Var(i2).Pos~=si)); continue; end
                Names{i2}=Var(i2).Name;
                ipl(i2,1)=i2;
                ipl(i2,2)=Var(i2).Plotno;
            end
        end
    end
    
    i2=find(ipl(:,1));
    ipl=ipl(i2,:);

    return;
end

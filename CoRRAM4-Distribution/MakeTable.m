function MakeTable(Mod, varargin)
%Produce an ASCII table of second moments from the respective arrays of Mod from simulation results

%{
    Copyright Alfred Mauﬂner

    Purpose: Produce a table of second moments from the solution and simulation
             of a DSGE model and write this table to an ASCII file. 

    Input (required):   Mod, and instance of the DSGE class.
    Input (optional):  Names, cell array of variable names
                       ip,    integer vector: the column indices of those variables for
                              which standard deviations (and first-order autocorrelations)
                              will be printed
                       ic,    integer vector: the colum indices of those variable for
                              which correlation coefficients will be printed.

    Remarks: the indices in ic must be a subset of the indices in ip.

    Revision history: 
        
         16 August    2019, first verion
         26 September 2019, error message added if correlation with variable not in the index set is required
         09 December  2020, message wrt employed policy function extended to allow for non-perturbation
                            solutions
         16 February  2022, flag Mod.Flags.Append checked, file will be appended without dialog, if
                            this flag is set.
         25 February  2023, check if number of variable names equals the number of print indices

   %}

% check existence of Corr_m
if isempty(Mod.Corr_m)
    error('The model has not been simulated yet.');
end

% check input arguments
if nargin==1 % gather information from Var structure
    [Table1,Table2,PNames,ColDes]=MakeTable1(Mod);
else
  PNames=varargin{1};
  ip=varargin{2};
   if max(size(PNames))<size(ip,1)
      error('The number of variable names must be equal to the number of print indices');
   end 
  ic=varargin{3};
  [Table1,Table2,ColDes]=MakeTable2(Mod,PNames,ip,ic);
end

% prepare print
fmtstr1='%6.2f';
fstring2=repmat(fmtstr1,1,size(Table1,2));
fstring=strcat('%',num2str(max(strlength(PNames))),'s',fstring2,'\n');

header=MakeHeader(Mod);

% open ASCII file: existing file will not be overwritten!
file=strcat(Mod.Outfile,'.txt');
% check if file exists
rc=exist(file,'file');
if rc==2
    if Mod.Flags.Append
        id=fopen(file,'a');
    else
        ok=questdlg(sprintf('File %s already exists. Overwrite?',file));
        if strcmp(ok,'Yes') || strcmp(ok,'Cancel')
            id=fopen(file,'w');
        else
            id=fopen(file,'a');
        end
    end
else
    id=fopen(file,'a');
end
fprintf(id,'%s\n','');
fprintf(id,'%s\n',datestr(datetime('now')));
fprintf(id,'%s\n','');
for ii=2:6
    fprintf(id,'%s\n',header{ii});
end
fprintf(id,'%s\n','');

np=length(PNames);

if ~Mod.Flags.Confidence
    for j1=1:np
         fprintf(id,fstring,PNames{j1},Table1(j1,:));
    end
else    
    fmstr1='    %6.2f      ';
    fmstr1=repmat(fmstr1,1,size(Table1,1));
    fmstr1=strcat('%',num2str(max(strlength(PNames))),'s',fmstr1);
    fmstr2=repmat('(%6.2f;%6.2f) ',1,size(Table1,2));
    fmstr2=strcat('%',num2str(max(strlength(PNames))),'s',fmstr2);
    for ii=1:np
        str1=sprintf(fmstr1,PNames{ii},Table1(ii,:));    
        str2=sprintf(fmstr2,' ',Table2(ii,:));
        fprintf(id,'%s\n',str1);
        fprintf(id,'%s\n',str2);
    end    
end
ntab=length(ColDes);
fprintf(id,'%s\n','');
for j1=1:1:ntab
    fprintf(id,'%s %s\n',ColDes{j1,:});
end

fclose(id);

if Mod.Flags.Excel
    file=strcat(Mod.Outfile,'.xlsx');
    [nr, nc]=size(Table1);
    Outc=cell(2*nr,1+nc);
    for ii=1:nr
        for jj=1:nc
            Outc{1+2*(ii-1),1}=PNames{ii};
            Outc{2+2*(ii-1),1}='';
            Outc{1+2*(ii-1),1+jj}=num2str(Table1(ii,jj),'%5.3f');
            Outc{2+2*(ii-1),1+jj}=strcat(num2str(Table2(ii,1+2*(jj-1)),'%5.3f'),'-',num2str(Table2(ii,2+2*(jj-1)),'%5.3f'));
        end
    end
    writecell(header,file);
    writecell(Outc,file,'Range','A8');
    writecell(ColDes,file,'Range',strcat('A',num2str(10+2*nr)));
end

return;

end

function [Table1,Table2,PNames,ColDes]=MakeTable1(Mod) % construct table from information provided by the Var structure
        
        
    % number of variables in the model
    nvar=Mod.nx+Mod.ny+Mod.nz;


    % initialize
    np=0;      % number of variables for which output is requested
    nr=0;      % number of variables for which relative standard deviations are requested
    nc=0;      % number of variables for which correlations are requested
    ind1=zeros(nvar,1); % indices for print
    ind2=zeros(nvar,1); % indices for relative deviations
    ind3=zeros(nvar,1); % indices for correlations
    PNames=cell(nvar,1);
    RNames=cell(nvar,1);
    CNames=cell(nvar,1);
        
    % gather information from Var
    for ii=1:1:nvar
        if Mod.Var(ii).Print>0
            np=np+1;        
            PNames{np}=Mod.Var(ii).Name;
            ind1(np)=ii;
            if Mod.Var(ii).Rel
                nr=nr+1;
                RNames{nr}=Mod.Var(ii).Name;
                ind2(nr)=ii;
             end
             if Mod.Var(ii).Corr
                nc=nc+1;
                CNames{nc}=Mod.Var(ii).Name;
                ind3(nc)=ii;
             end
        end
    end
    PNames=PNames(1:np);

    % extract information from Mod.Corr_m
    ColDes=cell(3+nr+nc,2); % cell with information about each column of the table
    ColDes{1,1}='Column  1 :';
    ColDes{1,2}='Variable';
    ColDes{2,1}='Column  2 :';
    ColDes{2,2}='Standard Deviation';
    string1='Column ';
    string2='Standard Deviation relative to Standard Deviation of Variable ';
    string3='Correlation with Variable ';

    tmp_m=Mod.Corr_m(:,:,1);
    tmp1=diag(tmp_m);
    
    % for further use: replace the main diagonal with ones to give the correct answer
    tmp_m=tmp_m-diag(diag(tmp_m))+eye(nvar);
    Table1=zeros(length(ind1(1:np)),2+nc+nr);
    if Mod.Flags.Confidence
        Table2=zeros(length(ind1(1:np)),2*(2+nc+nr));
    end

    Table1(:,1)=100*tmp1(ind1(1:np)); % standard deviations
    if Mod.Flags.Confidence
        tmp1=diag(Mod.Corr_l(:,:,1));
        Table2(:,1)=100*tmp1(ind1(1:np));
        tmp1=diag(Mod.Corr_u(:,:,1));
        Table2(:,2)=100*tmp1(ind1(1:np));
        tmp_l=Mod.Corr_l(:,:,1)-diag(diag(Mod.Corr_l(:,:,1)))+eye(nvar);
        tmp_u=Mod.Corr_u(:,:,1)-diag(diag(Mod.Corr_u(:,:,1)))+eye(nvar);
    else
        Table2=[];
    end
    for ii=1:nr
        Table1(:,1+ii)=Table1(:,1)./Table1(ii,1); % relative standard deviations
        if Mod.Flags.Confidence
            Table2(:,3+2*(ii-1))=Table2(:,1)./Table2(ii,1);
            Table2(:,4+2*(ii-1))=Table2(:,2)./Table2(ii,2);
        end
        ColDes{2+ii,1}=sprintf('%s %s %s',string1,num2str(2+ii,'%2i'),':');
        ColDes{2+ii,2}=sprintf('%s %s',string2,RNames{ii});     
    end
    for ii=1:nc
        Table1(:,1+nr+ii)=tmp_m(ind3(ii),ind1(1:np),1)'; % correlations
        if Mod.Flags.Confidence
            Table2(:,2*(1+nr)+1+(ii-1)*2)=tmp_l(ind3(ii),ind1(1:np),1)';
            Table2(:,2*(1+nr)+2+(ii-1)*2)=tmp_u(ind3(ii),ind1(1:np),1)';
        end
        ColDes{2+nr+ii,1}=sprintf('%s %s %s',string1,num2str(2+nr+ii),':');
        ColDes{2+nr+ii,2}=sprintf('%s %s',string3,CNames{ii});
    end

    tmp1=diag(Mod.Corr_m(:,:,2));
    Table1(:,2+nr+nc)=tmp1(ind1(1:np));
    if Mod.Flags.Confidence
        tmp1=diag(Mod.Corr_l(:,:,2));
        Table2(:,2*(1+nr+nc)+1)=tmp1(ind1(1:np));
        tmp1=diag(Mod.Corr_u(:,:,2));
        Table2(:,2*(1+nr+nc)+2)=tmp1(ind1(1:np));
    end
    
    ColDes{3+nr+nc,1}=sprintf('%s %s %s',string1,num2str(3+nr+ii),':');
    ColDes{3+nr+nc,2}='First Order Autocorrelation';
       
    return
end

function [Table1,Table2,ColDes]=MakeTable2(Mod,PNames,ip,ic)

    % prepare output array
    np=length(ip);
    nc=length(ic);
    ncol=nc+3;
    Table1=zeros(np,nc+2);
    
    ColDes=cell(ncol,2);
    for ii=1:ncol
        ColDes{ii,1}=sprintf('Column %i :',ii);
    end
    ColDes{1,2}='Variable name';
    ColDes{2,2}='Standard deviation';
    for ii=1:nc
        k=find(ip==ic(ii));
        if isempty(k); error('Variable no %i is not an element of print index',ic(ii)); return; end %#ok<UNRCH>
        ColDes{2+ii,2}=sprintf('Correlation with variable %s',PNames{k}); 
    end
    ColDes{ncol,2}='First order autocorrelation';

    % gather data
    tmp_m=Mod.Corr_m(:,:,1);
    nvar=size(tmp_m,2);

    tmp1=diag(Mod.Corr_m(:,:,1));
    Table1(:,1)=100*tmp1(ip);
    tmp_m=tmp_m-diag(diag(tmp_m))+eye(nvar);
    Table1(:,2:nc+1)=tmp_m(ip,ic,1);
    tmp1=diag(Mod.Corr_m(:,:,2));
    Table1(:,nc+2)=tmp1(ip);
    
    if Mod.Flags.Confidence
        Table2=zeros(np,2*(2+nc));
        tmp1=diag(Mod.Corr_l(:,:,1));
        Table2(:,1)=100*tmp1(ip);
        tmp1=diag(Mod.Corr_u(:,:,1));
        Table2(:,2)=100*tmp1(ip);
        tmp_l=Mod.Corr_l(:,:,1);
        tmp_l=tmp_l-diag(diag(tmp_l))+eye(nvar);
        tmp_u=Mod.Corr_u(:,:,1);
        tmp_u=tmp_u-diag(diag(tmp_u))+eye(nvar);
        ind1=3:2:(2*nc+2);
        ind2=4:2:(2*nc+2);        
        Table2(:,ind1)=tmp_l(ip,ic);
        Table2(:,ind2)=tmp_u(ip,ic);   
        tmp1=diag(Mod.Corr_l(:,:,2));
        Table2(:,2*nc+3)=tmp1(ip);
        tmp1=diag(Mod.Corr_u(:,:,2));
        Table2(:,2*nc+4)=tmp1(ip);
    else
        Table2=[];
    end

return

end

function header=MakeHeader(Mod)

    header=cell(6,1);
    header{1}=datestr(datetime('now'));
    if Mod.Flags.CheckBounds
        header{2}=sprintf('Second moments from %5i valid simulations of length = %5i and burnin = %5i.', Mod.valsim,Mod.nobs,Mod.burnin);
        header{3}='Bounds were checked';     
    else
        if Mod.Flags.Simulate
            header{2}=sprintf('Second moments from %5i simulations of length = %5i and burnin = %5i.',Mod.nofs,Mod.nobs,Mod.burnin);
            header{3}='Bounds were not checked';
        else
            header{2}='Analytical second moments.';
        end
    end

    header{4}=sprintf('Filter: %s',Mod.Filter);
    if strcmp('HP',Mod.Filter)
        header{5}=sprintf('Filter weight= %5.0f',Mod.hpl);
    else
        header{5}='';
    end
    
    if Mod.Flags.Simulate
        switch Mod.order
            case 1
                header{6}='Linear policy functions were used.';
            case 2
                if Mod.Prune==2
                    header{6}='Quadratic policy functions and the pruned system were used';
                else
                    header{6}='Quadratic policy functions were used.';
                end
            case 3
                header{6}='Cubic policy functions were used.';
            otherwise
                header{6}='Policy functions from other than perturbation solution were used';
        end
    end

    return
end
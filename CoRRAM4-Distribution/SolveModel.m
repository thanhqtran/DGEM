function M = SolveModel(M)
%  Computes the matrices of the first-, second-, and third-order perturbation solution of a DSGE model
%{
    Copyright: Alfred Maußner
    
    Purpose:   This function computes the matrices that make up the perturbation solution
               o the canonical DSGE model as defined in Heer and Maußner, Dynamic General
               Equilibrium Modeling, 3rd Ed., Chapter 2, Section 5.

    Revision history:

        13 March 2019, first version (from a previous version of SolveModel) 
        27 March 2019, extended to 3rd order part of solution
        24 May   2019, ordering of variables changed to (x(t+1),y(t+1),z(t+1),x(t),y(t),z(t))

    Input: M, an instance of the DSGE class

    Ouput: M, this instance with solution added to this instance

%}

% check existence of the model's equations
M.rc=M.ExistsModel();
if M.rc~=2; M.rc=1; return; end

% create and open logfile
M.lfh=fopen(strcat(M.Outfile,'_LogFile.txt'),'w');
fprintf(M.lfh,strcat(char(datetime('now')),'\n'));

% compute the derivative matrices from the model's equations
if strcmp(M.Derivatives,'ND')
    [M.FN, M.DN, M.HN,M.rc]=GetFDH_ND(M);
elseif strcmp(M.Derivatives,'AD')
    [M.FN, M.DN, M.HN, M.TN, M.rc]=GetFDHT_AD(M);
elseif strcmp(M.Derivatives,'SD')    
    [M.FN,M.DS,M.DN,M.HS,M.HN,M.TS,M.TN,M.rc]=GetFDHT_SD(M);
end

if M.rc>0; return; end

% solve for the linear part
[M.Hxw, M.Hyw, M.rc]=Linear(M);
   
if M.rc>0; return; end
 
if M.order>1
    [M.Hxww, M.Hxss, M.Hyww, M.Hyss, M.rc]=Quadratic(M);
    if M.rc>0; return; end
end

if M.order>2
    [M.Hxwww,M.Hywww,M.Hxssw,M.Hyssw,M.Hxsss,M.Hysss,M.rc]=Cubic(M);
end

return;

end
function cy = HPF(y,mu )
%Hodrick-Prescott filter

%{

    Copyright Alfred Mauﬂner
    First version: 20 January 2017 (developed from the Fortran code hpfilter)
    Revision:      25 August  2018 (extended to matrix of data)

    Purpose:  returns the cyclical component of a time series y
              obtained from applying the Hodrick-Prescott filter to this series.

    Usage: cy=HPFilter(y,mu);

    Input:   y : nobs by nvar, vector with observations

            mu : scalar, the filter weight (i.e. mu=1600 for quarterly data)

    Output   cy: nobs by nvar, vector, the cyclical components of y

    Remarks: This function call the Lapack routine DPBSV via the matlab function lapack.
             This allows me to just set up a 3 by size(y,1) matrix of coefficients to
             solve for trend instead of a full size(y,1) by size(y,1) matrix as does the
             function hpfilter from the econometrics toolbox.

%}

% determine size of y and check for proper input
[nobs nvar]=size(y);
if nobs<=3; error('Input must be an n by m matrix with n>3'); end;
if ~isreal(y); error('HP-Filter cannot process complex numbers'); end; 
% set up the upper 3 diagonals
m=zeros(3,nobs);

for i=1:nobs;
    m(1,i)=mu;
    m(2,i)=-4*mu;
    m(3,i)=1+6*mu;
end
   
m(1,1)=0;
m(1,2)=0;
m(2,1)=0;
m(2,2)=-2*mu;
m(2,nobs)=-2*mu;
m(3,1)=1+mu;
m(3,2)=1+5*mu;
m(3,nobs)=1+mu;
m(3,nobs-1)=1+5*mu;    

% parameters for calling DPBSV
UPLO='U';
KD=2;
NRHS=1;
LDAB=KD+1;
INFO=0;
cy=zeros(nobs,nvar);
for i=1:nvar;
    result=lapack('DPBSV',UPLO,nobs,KD,NRHS,m,LDAB,y(:,i),nobs,INFO);

    % check for correct solution
    INFO=result{9};
    if INFO ~=0; error('Not able to solve obtain the trend'); end;
    trend=result{7};
    cy(:,i)=y(:,i)-trend;
    
end;
  
return
end
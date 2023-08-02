function xf = HamiltonFilter(x,varargin)
% Computes the cyclical component of time series from the filter proposed by Hamilton (2017)

%{

    Alfred Maußner
    
    Purpose: computes the cyclical component of time series from the filter
             proposed by James Hamilton, Why You Should Never Use the Hodrick-Prescott Filter,
             NBER Working Paper No. W23429, 2017

    Input:   x, nobs by nvar matrix, each column stores observations on one of the nvar
                variables

             k, scalar, optional, the number of observations to discard, if left out, set to k=8
                as proposed by Hamilton for quarterly data.

    Output:  xf, nobs-k by nvar matrix, the cyclical components.


    Revision history:

        18 April     2019, first version from the procedure in BK_Filtered_Data.g
        29 September 2022, revised, regression is: y(t+h) on constant, y(t), y(t-1), y(t-2), y(t-3)!

%}

% get information
[nobs, nvar]=size(x);
if nargin>1
   h=varargin{1};
else
   h=8;
end
xf=zeros(nobs-h-3,nvar);

% loop over variables
for v=1:nvar
    xmat=[ones(nobs-h-3,1),x(4:nobs-h,v),x(3:nobs-h-1,v),x(2:nobs-h-2,v),x(1:nobs-h-3,v)];    
    b=xmat\x(h+4:nobs,v);
    xf(:,v)=x(h+4:nobs,v)-xmat*b;
end
return;
end
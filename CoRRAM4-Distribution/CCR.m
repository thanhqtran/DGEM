function CorMat = CCR(x,kbar)
%Computes an array of correlations between the columns of x for
%lags k=0, ..., kbar

%{
    Alfred Mauﬂner

    Purpose: compute a cross correlation matrix at lags k=0, 1, ..., kbar
	         from a data matrix x with nobs observations on nvar time series.

    Input:  x, nobs by nvar matrix
          kbar, scalar the maximim number of lags considered

    Output: CorMat, nvar by nvar by 1+kbar array, CorMat(i,j,k) stores
            the correlation coefficients between x(t,i) and x(t-k,j).
            The diagonal elements of CorMat(:,:,1) store the standard deviations.
 
    Remarks: the formula for the covariance is taken from the EVwies user guide and was tested
	         see the file CovTest.g

    Revision history:

        11 April 2019, first version, from the Gauss procedure CCR

			 
%}

% dimensions

[nobs, nvar]=size(x);
	
% initialize
CorMat=zeros(nvar,nvar,1+kbar);
xbar=mean(x);
    
% compute covariances
for k=1:1+kbar
    for ii=1:nvar
        for jj=1:nvar
                CorMat(ii,jj,k)=(x(k:nobs,ii)-xbar(ii))'*(x(1:nobs+1-k,jj)-xbar(jj))/nobs;
        end
    end
end
Sd=sqrt(diag(CorMat(:,:,1)));
D=diag(1./Sd);
for k=1:1+kbar
    CorMat(:,:,k)=D*CorMat(:,:,k)*D;
end
CorMat(:,:,1)=CorMat(:,:,1)+diag(Sd-diag(CorMat(:,:,1)));

return;
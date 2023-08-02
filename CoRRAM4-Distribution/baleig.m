function [Dl, Dr, iter] = baleig(A,B,max_iter)
% Performs two-sided scaling Dl\A*Dr, Dl\B*Dr in order to improve
% the sensitivity of generalized eigenvalues.

%{ 
    Alfred Maußner

    20 August 2019
    
    Code and description taken from Lemonnier and Van Dooren, Balancing
    regular matrix pencils, SIAM J. Matrix Anal. Appl., Vol. 28, 2006, pp. 253-263, p. 662.

    The diagonal matrices Dl and Dr are constrained to powers of 2 and are computed iteratively
    until the number of iterations max_iter is met or until the norms are
    between 1/2 and 2. Convergence is often reached after 2 or 3 steps.
    The diagonals of the scaling matrices are returned in Dl and Dr
    and so is iter, the number of iterations steps used by the method.

%}
n=size(A,1);
Dl=ones(1,n);
Dr=ones(1,n);
M=abs(A).^2+abs(B).^2;
for iter=1:max_iter,emax=0;emin=0;
    for i=1:n
        % scale the rows of M to have approximate row sum 1
        d=sum(M(i,:));e=-round(log2(abs(d))/2);
        M(i,:)=pow2(M(i,:),2*e);
        % apply the square root scaling also to Dl
        Dl(i)=pow2(Dl(i),-e);
        if e > emax, emax=e; end; if e < emin, emin=e; end
    end
    for i=1:n
        % scale the columns of M to have approximate column sum 1
        d=sum(M(:,i));e=-round(log2(abs(d))/2);
        M(:,i)=pow2(M(:,i),2*e);
        % apply the square root scaling also to Dr
        Dr(i)=pow2(Dr(i),e);
        if e > emax, emax=e; end; if e < emin, emin=e; end
    end
    % Stop if norms are all between 1/2 and 2
    if (emax<=emin+2), break; end
end

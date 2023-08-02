function Jac = CDJacRE(f,x0 )
% Computes the Jacobian matrix using central differences and Richardson extrapolation
% Usage: Jac=CDJaREc(f,x0)

%{ 
    Copyright: Alfred Maußner

    Purpose: computes the Jacobian matrix, i.e., the matrix of first-order partial derivatives
    of a user supplied, vector-valued function f:R^n->R^m at the point x0 using
    the central difference formula and Richardson extrapolation.

    Input: f, pointer to the function f. The function must return the m by 1 vector [f1; f2; ...;fm]=f(x0)
          x0, n by 1 vector, point at which f is to be evaluated

    Output: Jac, m by n matrix, where the element Jac(i,j) stores the partial derivative of the
            i-th element function fi with respect to variable x0(j).

    Remarks: See Herr and Maußner, 3rd. ed. for the respective algorithm.

    Revision history: 01 March 2019, first version

          
%}   

% parameters of the algorithm
h=0.5;
kmax=10;
eps1=eps;

% get the function value to determine m
f0=f(x0);

n=length(x0);
m=length(f0);

Jac=zeros(m,n);

vh=h*2.^(-(0:kmax));
f1=zeros(m,kmax+1);

for i=1:n; % loop over the elements of x0
    dvec=max([x0(i) h]).*vh;
    x1=repmat(x0,1,kmax+1);
    x2=x1;
    x1(i,:)=x1(i,:)+dvec;
    x2(i,:)=x2(i,:)-dvec;
    for j=1:kmax+1;
        f1(:,j)=(f(x1(:,j))-f(x2(:,j)))/(2*dvec(j)); % central differences
    end;
    for j=1:m;
        if sum(f1(j,:)~=0)>0; Jac(j,i)=Der(f1(j,:)',kmax,eps1); end;
    end;
end;

return;

    function est=Der(v,k,eps1) % implements Algorithm 13.3.1
        D=zeros(k+1,2);
        D(:,1)=v;
        err1=1;
        l=1;
        while l<=k;
            D(1:k+1-l,2)=D(2:k+2-l,1)+(D(2:k+2-l,1)-D(1:k+1-l,1))/(4^l-1);
            err2=2*abs(D(k+1-l,2)-D(k+2-l,1))/(abs(D(k+1-l,2))+abs(D(k+2-l,1))+eps1);
            if err2<err1; err1=err2; D(:,1)=D(:,2); else break; end;
            l=l+1;
        end;
        if l<k; est=D(k+2-l,1); else est=D(1,2); end;
        
        return;
    end
    
end


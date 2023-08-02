function Hmat = CDHesse(f,x0)
% Computes the Hesse matrix using the central difference and a four-point formula.
% Usage: Hmat=CDHesse(f,x0)

%{ 
    Copyright: Alfred Mauﬂner

    Purpose: computes the Hessian matrix, i.e., the matrix of second-order partial derivatives
    of a user supplied, vector-valued function f:R^n->R at the point x0 using
    the central difference formula for the second-order derivatives and a four-point
    formula for the mixed second-order derivatives.

    Input: f, pointer to the function f.
          x0, n by 1 vector, point at which f is to be evaluated

    Output: Hmat, n by n matrix, where the element Hmat(i,j) stores the second-order partial derivative of 
            f with respect to x0(i) and x0(j)

    Remarks: The diagonal elements of Hmat are computed from the central difference formula,
             the off diagonal elements are computed from a four point formula. 
             See Heer and Mauﬂner, 3rd. ed., Section 13.3.1.

    Revision history: 01 March 2019, first version
                      15 March 2019, since sign(0)=0, check for x0(i)==0 added
          
%}   

n=length(x0);
Hmat=zeros(n,n);
h=zeros(n,1);
eps1=eps^(1/3);
f0=f(x0); % function value at x0, used for central difference formula

%step sizes
for i=1:n;
    if x0(i)~=0; h(i)=sign(x0(i))*eps1*max(abs([x0(i) 1])); else h(i)=eps1; end;
end;
temp=x0+h;
h=temp-x0; % see Dennis and Schnabel (1983), p. 321 for this trick
ee=diag(h);

for i=1:n; % outer loop over x0(i)
    Hmat(i,i)=(f(x0+ee(:,i))+f(x0-ee(:,i))-2*f0)/(h(i)^2);
    for j=i+1:n; % innter loop over x0(j)
        Hmat(i,j)=(f(x0+ee(:,i)+ee(:,j))+f(x0-ee(:,i)-ee(:,j))-f(x0+ee(:,i)-ee(:,j))-f(x0-ee(:,i)+ee(:,j)))/(4*h(i)*h(j));
        Hmat(j,i)=Hmat(i,j);
    end;
end;
return;
end


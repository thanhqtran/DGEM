function Jac = CDJac(f,x0)
% Computes the Jacobian matrix using the central difference formula.
% Usage: Jac=CDJac(f,x0)

%{ 
    Copyright: Alfred Mauﬂner

    Purpose: computes the Jacobian matrix, i.e., the matrix of first-order partial derivatives
    of a user supplied, vector-valued function f:R^n->R^m at the point x0 using
    the central difference formula.

    Input: f, pointer to the function f. The function must return the m by 1 vector [f1; f2; ...;fm]=f(x0)
          x0, n by 1 vector, point at which f is to be evaluated

    Output: Jac, m by n matrix, where the element Jac(i,j) stores the partial derivative of the
            i-th element function fi with respect to variable x0(j).

    Remarks: Jac(i,j) is computed from the central difference formula. See Heer and Mauﬂner, 3rd. ed.
             Section 13.3.1.

    Revision history: 01 March 2019, first version

          
%}   

f0=f(x0); % to be able to determine the size of Jac
n=length(x0);
m=length(f0);

Jac=zeros(m,n);

% step size
eps1=eps^(1/3);
x1=x0;
x2=x0;

% loop over the elements of x0
for i=1:n;
    if x0(i)~=0; h=sign(x0(i))*eps1*max(abs([x0(i) 1])); else h=eps1; end;
    temp=x0(i);
    x1(i)=temp+h;
    x2(i)=temp-h;
    h=x1(i)-temp; % trick, see Dennis and Schnabel (1983), p. 99
    Jac(1:m,i)=(f(x1)-f(x2))/(2*h);
    x1(i)=x0(i);
    x2(i)=x0(i);
end;

return;

end


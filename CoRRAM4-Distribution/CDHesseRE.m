function Hmat = CDHesseRE(f,x0)
% Computes the Hesse matrix using the central difference and a four-point formula and Richardson extrapolation.
% Usage: Hmat=CDHesseRE(f,x0)

%{ 
    Copyright: Alfred Maußner

    Purpose: computes the Hessian matrix, i.e., the matrix of second-order partial derivatives
    of a user supplied, vector-valued function f:R^n->R at the point x0 using
    the central difference formula for the second-order derivatives and a four-point
    formula for the mixed second-order derivatives.

    Input: f, pointer to the function f.
          x0, n by 1 vector, point at which f is to be evaluated

    Output: Hmat, n by n matrix, where the element Hmat(i,j) stores the second-order partial derivative of 
            f with respect to x0(i) and x0(j)

    Remarks: See Herr and Maußner, 3rd. ed. for the respective algorithm.

    Revision history: 01 March 2019, first version

          
%}   

% parameters of the algorithm
h=0.5;
kmax=9;
eps1=eps;

f0=f(x0); % function value at x0, used for central difference formula

n=length(x0);

Hmat=zeros(n,n); % initialize Hmat

vh=h*2.^(-(0:kmax)); % steps

f1=zeros(kmax+1,1); % required for extrapolation

for i=1:n; % outer loop over elements of x0(i)
    dvec1=max([x0(i) h]).*vh;
    x11=repmat(x0,1,kmax+1);
    x22=x11;
    x11(i,:)=x11(i,:)+dvec1;
    x22(i,:)=x22(i,:)-dvec1;
    for j=1:(kmax+1);
        f1(j)=(f(x11(:,j))+f(x22(:,j))-2*f0)/(dvec1(j)^2);
    end;
    if sum(f1~=0)>0; Hmat(i,i)=Der(f1,kmax,eps1,2); end;
        
    for k=i+1:n; % inner loop over elements of x0(k)
        dvec2=max(abs([x0(k) h])).*vh;
        x1=x11;
        x2=x22;
        x3=x22;
        x4=x11;
        x1(k,:)=x1(k,:)+dvec2; % x(i)+h(i),x(j)+h(j)
        x2(k,:)=x2(k,:)-dvec2; % x(i)-h(i),x(j)-h(j)
        x3(k,:)=x3(k,:)+dvec2; % x(i)-h(i),x(j)+h(j)
        x4(k,:)=x4(k,:)-dvec2; % x(i)+h(i),x(j)-h(j)
        for j=1:(kmax+1);
            f1(j)=(f(x1(:,j))+f(x2(:,j))-f(x3(:,j))-f(x4(:,j)))/(4*dvec1(j)*dvec2(j));
        end;
        if sum(f1~=0)>0; Hmat(i,k)=Der(f1,kmax,eps1,4); Hmat(k,i)=Hmat(i,k); end;
    end;
end;
return;


 function est=Der(v,k,eps1,s) % implements Algorithm 13.3.1 with s={2,4}
        D=zeros(k+1,2);
        D(:,1)=v;
        err1=1;
        l=1;
        while l<=k;
            D(1:k+1-l,2)=D(2:k+2-l,1)+(D(2:k+2-l,1)-D(1:k+1-l,1))/(s^l-1);
            err2=2*abs(D(k+1-l,2)-D(k+2-l,1))/(abs(D(k+1-l,2))+abs(D(k+2-l,1))+eps1);
            if err2<err1; err1=err2; D(:,1)=D(:,2); else break; end;
            l=l+1;
        end;
        if l<k; est=D(k+2-l,1); else est=D(1,2); end;
        
        return;
    end
   

end


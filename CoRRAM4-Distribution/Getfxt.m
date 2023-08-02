function fx = Getfxt( x )
%Computes the matrix fxtilde from x
%{
    Author: Alfred Mauﬂner

    Date:   27 March 2019

    Purpose: Given the Jacobi matrix x of the vector-valued function
             f(x) with component functios fi(x), i=1, ..., m,
             this function computes the matrix

             fx:= [kron(eye(n),x(1,:); ..., kron(ey(n),x(m,:)],

             where n is the number of arguments of fi(x).
 %}

[m,n]=size(x);
fx=zeros(m*n,n^2);
for ii=1:m;
    fx(1+(ii-1)*n:ii*n,:)=kron(eye(n),x(ii,:));
end;
return;
end


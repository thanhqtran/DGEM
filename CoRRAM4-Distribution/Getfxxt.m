function fx = Getfxxt( x )
%Computes the matrix fxxtilde from x
%{
    Author: Alfred Mauﬂner

    Purpose: Given a matrix of nb blocks of size n by n,
             this function returns the nb times n^2 matrix
             whose ith row is the ith block arranged in one row, i.e.,

             fx(i,:)=vec(x(1+(i-1)*n:i*n,:)'

    Revision history:
 
    27 March 2019, first version
    29 March 2019, block size as optional input added

 %}
[m,n]=size(x);
nb=m/n;

fx=zeros(nb,n^2);
for ii=1:nb;
    fx(ii,:)=reshape(x(1+(ii-1)*n:ii*n,:),1,n^2);
end;
return;
end


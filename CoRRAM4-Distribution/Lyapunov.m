function X = Lyapunov(Amat,Bmat)
%Solve Lyapunov equation 
%{
    Alfred Mauﬂner
    22 February 2017, first version

    This program solves the Lyapunov equation

    X = Amat*X*Amat^T + Bmat.

    For this purpose, the equation is written as generalized Sylvester equation

    A*R-L*B=C
    D*R-L*E=F

    where: A=eye(n)
           B=Amat^T
           C=Bmat
           D=Amat
           E=eye(n)
           F=zeros(n,n)

    The Sylvester equation is solved by calling my function Sylvester. This function
    invokes the LAPACK routine DTGSYL to solve the equation.

%}

n=size(Amat,1);

X=Sylvester(eye(n),Amat',Bmat,Amat,eye(n),zeros(n,n));

return;

end


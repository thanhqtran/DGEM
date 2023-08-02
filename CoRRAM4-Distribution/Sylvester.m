function R = Sylvester(A,B,C,D,E,F);
% Solve the generalized Sylvester equation

%{
    Copyright Alfred Mauﬂner

    First version: 11 July 2015

    Purpose: Solve the generalized Sylvester equation

    Remarks: The function calls the function lapack which provides a Maltab interface
             to the Lapack function DTGSYL

     AR - LB = C
     DR - LE = F
 
     where R and L are unknown and A, B, C, D, E and F are given matrices. The dimensions are:

     R and L are m by n
     A and D are m by m
     B and E are n by n
     C and F are m by n.
    
     The steps to solve this equation are:

     Step 1:
     The QZ-Factorization of the pencils (A,D) and (B,E) is computed via the Matlab command qz giving

     A1=Q1*A*Z1
     D1=Q1*D*Z1

     B1=Q2*B*Z2
     E1=Q2*E*Z2

     Step 2: A1, B1, Q1*C*Z2, D1, E1, and Q1*F*Z2 are passed to the Lapack routine DTGSYL which returns R1

     Step 3: R=Z2*R1*V2'.

%}

%     find n and m
[m,i]=size(A);
[n,i]=size(B);

% Step 1:
[A1,D1,Q1,Z1] = qz(A,D,'real');
[B1,E1,Q2,Z2] = qz(B,E,'real');
%
% Step 2:
C1=Q1*C*Z2;
F1=Q1*F*Z2;

trans='N';
ijob=0;
scale=1;
dif=0;
work=0;
lwork=1;
iwork=zeros(n+m+6);
info=0;
result=lapack('DTGSYL',trans,ijob,m,n,A1,m,B1,n,C1,m,D1,m,E,n,F1,m,scale,dif,work,lwork,iwork,info);
lwork=result{20};
info=result{22};
R=result{9};
R=Z1*R*transpose(Z2);

return;

end


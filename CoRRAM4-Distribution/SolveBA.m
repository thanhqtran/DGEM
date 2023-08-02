function [Hx_w,Hy_w,rc] = SolveBA(Amat,Bmat,nx,ny,nz,lfh,balance,i_tol) 
% Obtain linear part of perturbation solution of DSGE models via QZ decomposition

%{
    Copyright:  Alfred Maußner    

    Purpose:  obtain the coefficients of the linear solution of a DSGE model via the QZ decomposition

              The model is given by

              (M): Bmat*E_t([(w(t+1)-w);(y(t+1)-v)]=Amat*([(w(t)-w);(y(t)-y)];

              where w(t):=[x(t);z(t)]

              and its solution is given by

              x(t+1)=Hx_w*(w(t)-w)
              y(t)  =Hy_w*(w(t)-w)

              (M) is the linearized version of the canonical model
              
              (CM1): 0=E_t g(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))
              (CM2): 0=z(t+1)-Rho*z(t)+Omega*epsilon(t+1).

    Revision history:
          13 March   2019 (from previous versions of the algorithm)
          20 August  2019, balancing of the matrix pencil added

    Input:  Amat, nw by nw matrix
            Bmat, nw by nw matrix,
              nx, integer, number of elements of the vector x(t)
              ny, integer, number of elements of the vector y(t)
              nz, integer, number of elements of the vector z(t)
             lfh, file handle for the log-file
           i_tol, scalar real, tolerance in eliminating imaginary parts of the complex solution

    Output:  Hx_w, nx by (nx+nz) matrix, the linear part of the solution for the model's states
             Hy_w, ny by (nx+nz) matrix, the linear part of the solution for not predetermined variables
             rc,   integer, the index of the error messages, see the code in SolveModel for these messages

%}


% Initialize
warning('off','MATLAB:singularMatrix');
lastwarn('');
rc=0;
Hx_w=[];
Hy_w=[];
nvar=nx+ny+nz;

% obtain the complex, ordered generalized Schur decomposition of the pencil (Amat, Bmat)
if balance
     [dl,dr,~]=baleig(Amat,Bmat,10);
    % ab=zeros(nvar,nvar);
    % bb=zeros(nvar,nvar);
    % for i=1:nvar
    %    for j=1:nvar
    %        ab(i,j)=(dr(j)/dl(i))*Amat(i,j);
    %        bb(i,j)=(dr(j)/dl(i))*Bmat(i,j);
    %    end
    % end
    ab=diag(1./dl)*Amat*diag(dr);
    bb=diag(1./dl)*Bmat*diag(dr);
    % ok1=max(max(abs(abt-ab)));
    % ok2=max(max(abs(bbt-bb)));
    [T,S,Q,Z]=qz(ab,bb,'complex');
    [T,S,~,Z]=ordqz(T,S,Q,Z,'udi');
else
    [T,S,Q,Z] = qz(Amat,Bmat,'complex');
    [T,S,~,Z] = ordqz(T,S,Q,Z,'udi');
end
lambda=ordeig(T,S);
out=[1:length(Amat);lambda'];
fprintf(lfh,'Eigenvalue no %3i = %10.6f\n',out);

% test for stability
test1=sum(abs(lambda(1:nx+nz))>1)>0;
if test1
    fprintf(lfh,'Not all of the %4i eigenvalues are within the unit circle\n',nx+nz);
    rc=5;
    return;
end

% test for determinancy
test1=sum(abs(lambda(nx+nz+1:nx+nz+ny))<1)>0;
if test1
    fprintf(lfh,'Not all of the %4i eigenvalues are outside the unit circle\n',nv);
    rc=6;
    return;
end
Temp1=linsolve(S(1:nx+nz,1:nx+nz)',Z(1:nx+nz,1:nx+nz)');
if MatrixSingular; rc=4; return; end
Temp2=linsolve(Z(1:nx+nz,1:nx+nz)',T(1:nx+nz,1:nx+nz)');
if MatrixSingular; rc=4; return; end
Hw_w=(Temp1')*(Temp2');
Hy_w=linsolve(Z(1:nx+nz,1:nx+nz)',Z(nx+nz+1:nx+nz+ny,1:nx+nz)');
Hy_w=Hy_w';
test1=max(abs(imag(Hw_w)));
if test1>i_tol
    fprintf(lfh,'There are complex coefficients ind Hw_w: %10.6f\n',test1);
end
Hw_w=real(Hw_w);

test1=max(abs(imag(Hy_w)));
if test1>i_tol
    fprintf(lfh,'There are complex coefficients ind Hy_w: %10.6f\n',test1);
end
Hy_w=real(Hy_w);

if balance        
    Hw_w=diag(dr(1:nx+nz))*Hw_w*diag(1./dr(1:nx+nz));
    Hy_w=diag(dr(1+nx+nz:nvar))*Hy_w*diag(1./dr(1:nx+nz));
end

if nx>0; Hx_w=Hw_w(1:nx,:); end

function rc=MatrixSingular()
        rc=0;
        [~,msgid]=lastwarn;
        if ~isempty(msgid)
            if strcmp(msgid,'MATLAB:singularMatrix')
                rc=true;
            end
        end
    end
end


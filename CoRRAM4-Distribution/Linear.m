function [Hx_w,Hy_w,rc] = Linear(Mod) 
% Obtain linear part of perturbation solution of the canonical DSGE model

%{
    Copyright:  Alfred Maußner    

    Purpose:  obtain the coefficients of the linear solution of a DSGE model given by

              (CM1): 0=E_t g(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))
              (CM2): 0=z(t+1)-Rho*z(t)+sigma*Omega*epsilon(t+1).  

              Depending on the setting Mod.Flags.Reduce the linearized model is either given
              by 

              (BA): Bmat*E_t([(w(t+1)-w);(y(t+1)-v)]=Amat*([(w(t)-w);(y(t)-y)]

              or by

              (I): Cy*(y(t)-y) + Cx*(x(t)-x) + Cz*z(t)=0

              (II): Dx*E(x(t+1)-x) + Dz*E(z(t+1)) + Dy*E(y(t+1)-y) + Fx*(x(t)-x) + Fz*z(t) + Fy*(y(t)-y)=0
 
              This model is reduced in two steps. First, the matrix Cy is partitioned into Cu and Cv,
              where Cu consists of the first n(u) colums of Cy. This gives
 
                
              (I'): Cu*(u(t)-u) + Cv*(v(t)-v) + Cx*(x(t)-x) + Cz*z(t)=0

              (II'): Dx*E(x(t+1)-x) + Dz*E(z(t+1)) + Du*E(u(t+1)-y) +Dv*E(v(t+1)-v)
                     + Fx*(x(t)-x) + Fz*z(t) + Fu*(u(t)-u) + Fv*(v(t)-v)=0

              If Cu is singular a new matrix tCu and related matrices tCv, tDu, tDv, tFu, tFv
              are constructed from a QR decomposition with column pivoting of Cy:
 
              tCy:=Q^T*Cy*P=[tCu, tCv] giving 
  
              (I''): Q^T*Cy*P*(P^T y(t)-y) + Q^T Cx*(x(t)-x) + Q^t*Cz*z(t)=0
                    [tCu, tCv]*(u(t)-u|v(t)-v)) + Q^T Cx*(x(t)-x) + Q^t*Cz*z(t)=0

              where u(t) are now equal to the first n(u) rows of P^T(y(t)-y) and
              v(t) are equal to the remaining n(y)-n(u) rows.

              The columns of Dy and Fy are interchanged accordingly so that II reads
 
              (II'): Dx*E(x(t+1)-x) + Dz*E(z(t+1)) + tDu*E(u(t+1)-u) tDv*E(v(t+1)-v)
                     + Fx*(x(t)-x) + Fz*z(t) + Fu*(u(t)-y) + Fv*(v(t)-v)=0
              equations of the model.

              In the second step, (I') or (I'') are solved for u(t) so that (II') or  (II'')
              can be reduced to the model

              B*E((x(t+1)-x)|z(t+1)|(v(t+1)-v))+A*((x(t)-x)|z(t)|(v(t)-v))

              (where the definitions of B and A integrate the shock process).

              If B is invertible so that W=B^(-1)*A, the solution is obtained via
              the Schur decompositition of W. Otherwise, the solution is obtained
              via the factorization of the pencil (B,A) via a call to the function SolveBA.
              The solution of the model are the matrices in

              x(t+1)=Hx_w*(x(t);z(t))
              y(t)  =Hy_w*(x(t);z(t))


    Revision history:
           12 March     2019, first version 
           31 May       2019, do not reduce models with either nx==0 or nz==0
            3 September 2019, print absolute value of eigenvalues

    Input:  Model, an instance of the DSGE class

    Output:  Hx_w, nx by (nx+nz) matrix, the linear part of the solution for the model's states
             Hy_w, ny by (nx+nz) matrix, the linear part of the solution for not predetermined variables
             rc,   integer, the index of the error messages, see the code in SolveModel for these messages

    Remarks:  the program assumes the canonical model Eg(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))=0
              z(t+1)=Rho*z(t)+Omega*epsilon(t+1)

              The linearized model
 
               
%}

% Initialize 
Hx_w=[];
Hy_w=[];
rc=0;
P=[];

% and (for simplicity)
nx=Mod.nx;
nz=Mod.nz;
ny=Mod.ny;
nu=Mod.nu;
nv=ny-nu;
na=nx+nz+ny;

% models with either nx==0 or nz==0 will not be reduced
if nx==0 && nz==0
    error('nx and nz cannot be both equal to zero');
end
if nx==0 || nz==0
    if nx==0
        fprintf(Mod.lfh,'Hint: the model has no endogenous state and will not be reduced.');
    else
        fprintf(Mod.lfh,'Hint: deterministic models will not be reduced.');
    end
    Mod.Flags.Reduce=0;
end

% the BA model
if ~Mod.Flags.Reduce
    if nx>0
        if nz>0
            Amat=[-Mod.DN(:,1+na:2*na);[zeros(nz,nx), Mod.Rho, zeros(nz,ny)]];
            Bmat=[Mod.DN(:,1:na);[zeros(nz,nx),eye(nz),zeros(nz,ny)]];
        else
            Amat=-Mod.DN(:,1+na:2*na);
            Bmat= Mod.DN(:,1:na);    
        end
    else
        Amat=[-Mod.DN(:,1+na:2*na);[Mod.Rho, zeros(nz,ny)]];
        Bmat=[Mod.DN(:,1:na);[eye(nz),zeros(nz,ny)]];
    end
    [Hx_w, Hy_w, rc]=SolveBA(Amat,Bmat,nx,ny,nz,Mod.lfh,Mod.Flags.Balance,Mod.itol);
    return;
end

% the reduced model
% check for invertability of Cu and set up matrices
Cy=Mod.DN(1:nu,1+na+nx+nz:2*na);
if rank(Cy(:,1:nu))<nu
    fprintf(Mod.lfh,'Cu matrix is singular. The program proceeds with QR factorization.\n');
    [Q,R,P]=qr(Cy);
    Cu=R(:,1:nu);    
    Cv=R(:,1+nu:ny);
    Cx=Q'*Mod.DN(1:nu,1+na:na+nx);
    Cz=Q'*Mod.DN(1:nu,1+na+nx:na+nx+nz);
    Dy=Mod.DN(1+nu:nx+ny,1+nx+nz:na)*P;
    Du=Dy(:,1:nu);
    Dv=Dy(:,1+nu:ny);
    Fy=Mod.DN(1+nu:nx+ny,1+na+nx+nz:2*na)*P;
    Fu=Fy(:,1:nu);
    Fv=Fy(:,1+nu:ny);
else
    Cu=Cy(:,1:nu);
    Cv=Cy(:,1+nu:ny);
    Cx=Mod.DN(1:nu,1+na:na+nx);
    Cz=Mod.DN(1:nu,1+na+nx:na+nx+nz);
    Du=Mod.DN(1+nu:nx+ny,1+nx+nz:nx+nz+nu);
    Dv=Mod.DN(1+nu:nx+ny,1+nx+nz+nu:na);
    Fu=Mod.DN(1+nu:nx+ny,1+na+nx+nz:na+nx+nz+nu);
    Fv=Mod.DN(1+nu:nx+ny,1+na+nx+nz+nu:2*na);
end
Dx=Mod.DN(1+nu:nx+ny,1:nx);
Fx=Mod.DN(1+nu:nx+ny,1+na:na+nx);
Dz=Mod.DN(1+nu:nx+ny,1+nx:nx+nz);
Fz=Mod.DN(1+nu:nx+ny,1+na+nx:na+nx+nz);

% set up reduced system
CuiCwv=linsolve(Cu,[Cx Cz Cv]);
Bmat=[[Dx Dz Dv]-Du*CuiCwv;[zeros(nz,nx), eye(nz), zeros(nz,nv)]];
Amat=[-[Fx Fz Fv]+Fu*CuiCwv;[zeros(nz,nx), Mod.Rho, zeros(nz,nv)]];

if rank(Bmat)<(nx+nz+nv)
    fprintf(Mod.lfh,'Amat is singular\n');
    [Hx_w, Lvw, rc]=SolveBA(Amat,Bmat,nx,nv,nz,Mod.lfh,Mod.Flags.Balance,Mod.itol);
    if rc>0; return; end
else
    Wmat=linsolve(Bmat,Amat); 
    % obtain the complex, ordered Schur decomposition of Wmat
    [Tmat,Smat]=schur(Wmat,'complex');
    [Tmat,Smat]=ordschur(Tmat,Smat,'udi');
    lambda=ordeig(Smat);
    out=[1:length(Wmat);abs(lambda')];
    fprintf(Mod.lfh,'Eigenvalue no %3i = %10.6f\n',out);
    % test for stability
    test1=sum(abs(lambda(1:nx+nz))>1)>0;
    if test1
        fprintf(Mod.lfh,'Not all of the %4i eigenvalues are within the unit circle',nx+nz);
        rc=5;
        return;
    end
    % test for determinancy
    test1=sum(abs(lambda(nx+nz+1:nx+nz+nv))<1)>0;
    if test1
        fprintf(Mod.lfh,'Not all of the %4i eigenvalues are outside the unit circle',nv);
        rc=6;
        return;
    end
    Hw_w=linsolve(Tmat(1:nx+nz,1:nx+nz)',(Tmat(1:nx+nz,1:nx+nz)*Smat(1:nx+nz,1:nx+nz))');
    if MatrixSingular; rc=4; return; end
    Lvw=linsolve(Tmat(1:nx+nz,1:nx+nz)',Tmat(1+nx+nz:nx+nz+nv,1:nx+nz)');	
    Hw_w=Hw_w';
    Lvw=Lvw';
    test1=max(abs(imag(Hw_w)));
    if test1>Mod.itol
        fprintf(Mod.lfh,'There are complex coefficients ind Hw_w: %10.6f',itest);
    end
    Hw_w=real(Hw_w);
    test1=max(abs(imag(Lvw)));
    if test1>Mod.itol
        fprintf(Mod.lfh,'There are complex coefficients ind Lvw: %10.6f',itest);
    end
    Lvw=real(Lvw);    
    Hx_w=Hw_w(1:nx,:);
end

% final step compute remaining part of Hy_w
Luw=-CuiCwv*[eye(nx+nz);Lvw];
Hy_w=[Luw;Lvw];
if ~isempty(P)
    Hy_w=P*Hy_w;
end

return;

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


function [Table1,Table2,ColDes]=MakeTable2(Mod,PNames,ip,ic)

    % prepare output array
    np=length(ip);
    nc=length(ic);
    ncol=nc+3;
    Table1=zeros(np,nc+2);
    
    ColDes=cell(ncol,2);
    for ii=1:ncol
        ColDes{ii,1}=sprintf('Column %i :',ii);
    end
    ColDes{1,2}='Variable name';
    ColDes{2,2}='Standard deviation';
    for ii=1:nc
        k=find(ip==ic(ii));
        ColDes{2+ii,2}=sprintf('Correlation with variable %s',PNames{k}); %#ok<FNDSB>
    end
    ColDes{ncol,2}='First order autocorrelation';

    % gather data
    tmp_m=Mod.Corr_m(:,:,1);
    nvar=size(tmp_m,2);

    tmp1=diag(Mod.Corr_m(:,:,1));
    Table1(:,1)=100*tmp1(ip);
    tmp_m=tmp_m-diag(diag(tmp_m))+eye(nvar);
    Table1(:,2:nc+1)=tmp_m(ip,ic,1);
    tmp1=diag(Mod.Corr_m(:,:,2));
    Table1(:,nc+2)=tmp1(ip);
    
    if Mod.Flags.Confidence
        Table2=zeros(np,2*(2+nc));
        tmp1=diag(Mod.Corr_l(:,:,1));
        Table2(:,1)=100*tmp1(ip);
        tmp1=diag(Mod.Corr_u(:,:,1));
        Table2(:,2)=100*tmp1(ip);
        tmp_l=Mod.Corr_l(:,:,1);
        tmp_l=tmp_l-diag(diag(tmp_l))+eye(nvar);
        tmp_u=Mod.Corr_u(:,:,1);
        tmp_u=tmp_u-diag(diag(tmp_u))+eye(nvar);
        ind1=3:2:(2*nc+2);
        ind2=4:2:(2*nc+2);        
        Table2(:,ind1)=tmp_l(ip,ic);
        Table2(:,ind2)=tmp_u(ip,ic);   
        tmp1=diag(Mod.Corr_l(:,:,2));
        Table2(:,2*nc+3)=tmp1(ip);
        tmp1=diag(Mod.Corr_u(:,:,2));
        Table2(:,2*nc+4)=tmp1(ip);
    else
        Table2=[];
    end

return

end

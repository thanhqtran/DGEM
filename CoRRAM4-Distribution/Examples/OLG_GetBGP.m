function fx=OLG_GetBGP(x,Par)
%The zero of this function in x (for a given set of parameters Par) is the balanced growth path of t generation OLG model

[~,~,nvec,kvec,~] = OLG_BGP1(x,Par);

fx=ones(3,1);
fx(1)=sum(kvec)/Par.T-x(1);
fx(2)=sum(nvec)/Par.T-x(2);
fx(3)=kvec(1);

return;

end


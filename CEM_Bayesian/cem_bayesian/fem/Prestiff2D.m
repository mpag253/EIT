function [pK,XL_p_K,index]=Prestiff2D(Nodes,Elements,XL_indicator)
Ne=length(Elements);
ElementsT=Elements';
iix=repmat(reshape((1:9*Ne),9,Ne),3,1); iix=iix(:);
iiy=(ElementsT(:)*ones(1,9))';  iiy=iiy(:);
iiv=zeros(27,Ne);
indx=(ElementsT(:)*ones(1,3))';indx=indx(:);
indy=repmat(Elements',3,1); indy=indy(:);
L=[-1 1 0;-1 0 1];
for ii=1:Ne
    Jt=[(Nodes(Elements(ii,2),:)-Nodes(Elements(ii,1),:))',...
        (Nodes(Elements(ii,3),:)-Nodes(Elements(ii,1),:))']';
        N=Jt\L;
        Ints=1/6*abs(det(Jt))*(N'*N);
        iiv(:,ii)=[Ints(:);Ints(:);Ints(:)];
end
pK=sparse(iix,iiy,iiv(:));
index=[indx(:),indy(:)];

if XL_indicator
    Nn=length(Nodes);
    XL_p_K=pK(:);
    indF=find(XL_p_K~=0);
    mod_in=mod(indF,length(indy));
    mod_in(mod_in==0)=length(indy);
    Sindy=indy(mod_in);
    indxo=0:Nn-1; 
    indxo=ones(length(indx),1)*indxo;
    indxo=indxo(:);
    indxo=indxo(indF);
    indxo=Nn*indxo;
    Sindx=indx*ones(1,Nn);
    Sindx=Sindx(:);
    Sindx=Sindx(indF); 
    Sindx=Sindx+indxo;
    JJ=XL_p_K(indF);

    XL_p_K=sparse(Sindx,Sindy,JJ,...
        Nn^2,Nn,length(indF));
else
    XL_p_K=0;
end
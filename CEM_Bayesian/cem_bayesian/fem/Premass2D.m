function [pM,KL_pM,index]=Premass2D(Nodes,Elements,indicator)
Nn=length(Nodes);
Ne=length(Elements);
ElementsT=Elements';
iix=repmat(reshape((1:9*Ne),9,Ne),3,1); iix=iix(:);
iiy=(ElementsT(:)*ones(1,9))';  iiy=iiy(:);
indx=(ElementsT(:)*ones(1,3))';indx=indx(:);
indy=repmat(Elements',3,1); indy=indy(:);
iiv=zeros(27,Ne);
% I=[1/20 1/60 1/60; 1/60 1/60 1/120;1/60 1/120 1/60];
% I1=[1/20 1/60 1/60; 1/60 1/60 1/120;1/60 1/120 1/60];


% I1=1/12*[2 1 1;1 2 1;1 2 1];
I1=[1/20 1/60 1/60; 1/60 1/60 1/120;1/60 1/120 1/60];

Pe=[0 0 1; 1 0 0;0 1 0];
I2=Pe*I1*Pe';
I3=Pe*I2*Pe';
I1=I1(:);
I2=I2(:);
I3=I3(:);
J=zeros(1,Ne);
parfor ii=1:Ne
    Jt=[(Nodes(Elements(ii,2),:)-Nodes(Elements(ii,1),:))',...
        (Nodes(Elements(ii,3),:)-Nodes(Elements(ii,1),:))']';
%     Ints1=I*abs(det(Jt));    
%     Ints2=circshift(circshift(Ints1,1,1),1,2);
%     Ints3=circshift(circshift(Ints2,1,1),1,2);
%     iiv(:,ii)=[Ints1(:);Ints2(:);Ints3(:)];
J(ii)=abs(det(Jt));
end
iiv=kron([I1;I2;I3],J);
pM=sparse(iix,iiy,iiv(:));
index=[indx(:),indy(:)];

if indicator
    KL_pM=pM(:);
    indF=find(KL_pM~=0);
    Sindy=repmat(indy,1,Nn); 
    Sindy=Sindy(indF);
    indxo=0:Nn-1; 
    indxo=repmat(indxo,9*Ne,1);
    indxo=Nn*indxo(indF);
    Sindx=repmat(indx,1,Nn);
    Sindx=Sindx(indF)+indxo;
    KL_pM=sparse(Sindx,Sindy,KL_pM(indF),Nn^2,Nn,length(indF));
else
    KL_pM=0;
end
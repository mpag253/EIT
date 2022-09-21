function [pK,SpK,IND]=PreMassCurved1D(Nodes,Bdy_indicator,indicator)
NN=length(Nodes);
Node_indx=(1:NN)';

if Bdy_indicator
    Elements=[Node_indx(1:end-1),Node_indx(2:end);Node_indx(end) Node_indx(1) ];
else
    Elements=[Node_indx(1:end-1),Node_indx(2:end)];
end

NE=length(Elements); 
indx=zeros(4*NE,1); indy=indx;
iix=zeros(1,8*NE);iiy=zeros(1,8*NE);iiv=zeros(1,8*NE);

I=[1/4, 1/12; 1/12, 1/12];
for ii=1:NE
    forx=[1 1]'*Elements(ii,:);
    fory=Elements(ii,:)'*[1 1];
    indx(1+(ii-1)*4:ii*4)=forx(:);
    indy(1+(ii-1)*4:ii*4)=fory(:);
    JE=norm(Nodes(Elements(ii,2),:)-Nodes(Elements(ii,1),:));
    Ints=JE*I;
    iix(1+(ii-1)*8:ii*8)=[1+(ii-1)*4:4*ii,1+(ii-1)*4:4*ii];
    YYY=ones(4,1)*Elements(ii,:);
    iiy(1+(ii-1)*8:ii*8)=YYY(:);
    iiv(1+(ii-1)*8:ii*8)=[Ints(:);flipud(Ints(:))];
end
pK=sparse(iix,iiy,iiv);
IND=[indx(:),indy(:)];
SpK=pK(:);
if indicator
    Sindy=indy*ones(1,NN); Sindy=Sindy(:);
    indxo=0:NN-1; indxo=ones(length(indx),1)*indxo; indxo=NN*indxo(:);
    Sindx=indx*ones(1,NN); Sindx=Sindx(:)+indxo;
    indF=find(SpK~=0);
    SpK=sparse(Sindx(indF),Sindy(indF),SpK(indF),NN^2,NN,length(indF));
else
    SpK=0;
end
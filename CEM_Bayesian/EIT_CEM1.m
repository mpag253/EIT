clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example EIT code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('distmesh')
rng(1234)
iter=1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter=20; %Max number of Gauss-Newton iterations
h=.02; %mesh refinement

rad=1; %radius of the circular domain

n_elec=16;

n_curr_pats=n_elec;

cont_imp=0.001; %contact impedance (this is uniform over the electrodes)


L=n_elec;
z=cont_imp;

N=-speye(n_elec);
N(:,1)=[];
N(1,:)=1;

% what type of currents do we apply?
% 
% II=zeros(L,1);
% II(1)=1;
% II(3)=-1;

II=toeplitz([1;-1;zeros(n_elec-2,1)],[1,zeros(1,n_curr_pats-2),-1]);
MeasPattern=II;

% II=1/sqrt(2)*CurrentPattern(:,1:16);

% return

% what type of measurements can we get?
% n_curr_pats=16;

theta_elec=linspace(-pi,pi,2*n_elec+1)';
theta_elec=theta_elec(1:end-1);
elec_pts=[rad*cos(theta_elec),rad*sin(theta_elec)];

sig_0=1;


nl=1;
beta=1;
bb=1*1e1;
aa=1*1e0;
cc=4*1e1;
mu_sig=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mesh_shape=@(p)dcircle(p,0,0,rad);
mesh_refine=@(p) min(1-10*dcircle(p,0,0,rad),2);

pfix=elec_pts;
[nodes,~]=distmesh2d(mesh_shape,mesh_refine,h,[-rad,-rad;rad,rad],pfix);
nn=length(nodes);

dt=delaunayTriangulation(nodes);
tris=dt.ConnectivityList;
stop


x=nodes(:,1);
y=nodes(:,2);

bdy_el=freeBoundary(dt); %elements on boundary

full_bdy=x.^2+y.^2>rad^2-1e-6; %boundary tolerance
bdy_indx=find(full_bdy);
nn_b=sum(full_bdy);
theta_bdy=atan2(y(bdy_indx),x(bdy_indx));

[~,ss_b]=sort(atan2(y(full_bdy),x(full_bdy)));


in_elec=zeros(nn_b,n_elec);
n_per_elec=zeros(n_elec,1);
for ii=1:n_elec
    in_elec(:,ii)=theta_bdy>=theta_elec(2*ii-1)-eps & ...
        theta_bdy<=theta_elec(2*ii)+eps;  
    n_per_elec(ii)=sum(in_elec(:,ii));
    
end
in_elec=logical(in_elec);
min(theta_bdy)
max(theta_bdy)
theta_elec




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_true=ones(nn,1)-gaussian([-.5,.3],.3,.9,nodes)-...
    gaussian([.5,.3],.3,.9,nodes)+gaussian([0,-.5],.3,.9,nodes);


% sigma_true=ones(nn,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pK,SpK,i_K]=Prestiff2D(nodes,tris,1);
% K=sparse(i_K(:,1),i_K(:,2),pK*sigma_true,nn,nn);
K=sparse(i_K(:,1),i_K(:,2),pK*exp(sigma_true),nn,nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M_bdy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_B=sparse([],[],0,nn,nn);
for ii=1:n_elec
[pM_b,~,i_Mb]=PreMassCurved1D(nodes(bdy_indx(in_elec(:,ii)),:),0,0);

M_b=sparse(i_Mb(:,1),i_Mb(:,2),pM_b*ones(n_per_elec(ii),1),...
    n_per_elec(ii),n_per_elec(ii));
M_B(bdy_indx(in_elec(:,ii)),bdy_indx(in_elec(:,ii)))=M_b;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=K+1/z*M_B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elec_nodes=unique(bdy_el);

elec_nodes=bdy_el(:,1);

Chat=sparse(nn,L);
for ii=1:length(elec_nodes)
    for jj=1:L
        inter=0;
        if sum(elec_nodes(ii)==bdy_indx(in_elec(:,jj)))==1
            
            [rr,cc]=find(elec_nodes(ii)==bdy_el);
            
            if sum(bdy_el(rr(1),2-mod(cc(1)+1,2))==bdy_indx(in_elec(:,jj)))
                                    
                Chat(elec_nodes(ii),jj)=Chat(elec_nodes(ii),jj)+1/2*norm(...
                    nodes(bdy_el(rr(1),2-mod(cc(1)+1,2)),:)-...
                    nodes(bdy_el(rr(1),cc(1)),:));
            end
            if sum(bdy_el(rr(2),2-mod(cc(2)+1,2))==bdy_indx(in_elec(:,jj)))
                Chat(elec_nodes(ii),jj)=Chat(elec_nodes(ii),jj)+1/2*norm(...
                    nodes(bdy_el(rr(2),2-mod(cc(2)+1,2)),:)-...
                    nodes(bdy_el(rr(2),cc(2)),:));
            end
        end

    end
end
C=-1/z*Chat*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
el_lengths=zeros(L,1);
for ii=1:L
    el_length=0;
    indx=bdy_indx(in_elec(:,ii));
    for jj=1:length(bdy_indx(in_elec(:,ii)))-1
        el_length=el_length+norm(nodes(indx(jj+1),:)-nodes(indx(jj),:));
    end
    el_lengths(ii)=el_length;
end
D=1/z*N'*diag(el_lengths)*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=[zeros(nn,n_curr_pats);N'*II];

D
size(B)
size(C)
size(D)

A=[B,C;C',D];
uu=A\F;
stop
Meas_rest=sparse((1:n_elec-1),(nn+1:nn+n_elec-1),1);


data=MeasPattern*N*Meas_rest*uu;

data_v=data(:);




%%%%%%%%%%%%%%%%%%%% Inversions
% 
% [pK,~,i_K]=Prestiff2D(nodes,tris,0);
% K=sparse(i_K(:,1),i_K(:,2),pK*sigma_true,nn,nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bimp=1/z*M_B;





sig_new=ones(nn,1);

Le=1e+4*speye(n_curr_pats*n_elec);

prior_corr_length=.1;

G_sig=1*exp(-1/prior_corr_length^2*((x-x').^2+(y-y').^2))+1e-6*eye(nn);
L_sig=chol(inv(G_sig));

mu_sig=1;


while iter<=max_iter
    
    K=sparse(i_K(:,1),i_K(:,2),pK*exp(sig_new),nn,nn);
    B=K+Bimp;
    
    A=[B,C;C',D];
    uu_new=A\F;    
    
    
    d_L=MeasPattern*N*(Meas_rest/A);
    S=sparse(n_curr_pats*n_elec,nn);
    for jj=1:n_curr_pats
        S(1+(jj-1)*n_elec:jj*n_elec,:)=-d_L*[reshape(SpK*uu_new(1:nn,jj),nn,nn);zeros(n_elec-1,nn)];
    end
    J_total=[Le*S*diag(exp(sig_new));L_sig];
%     J_total=[Le*S;L_sig];

    b_total=[Le*(reshape(MeasPattern*N*Meas_rest*uu_new,n_curr_pats*n_elec,1)-data_v);...
        L_sig*(sig_new-mu_sig*ones(nn,1))];
    dir_k=-J_total\b_total;
    %%%%%%%%%%%%%%%%%%%%%%%%%% Cost Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P=@(kappa) norm(Le*(reshape(MeasPattern*N*Meas_rest*([sparse(i_K(:,1),i_K(:,2),pK*...
        exp(sig_new+kappa*dir_k))+Bimp,C;C',D]\F),n_curr_pats*n_elec,1)-data_v))^2+...
        norm(L_sig*(sig_new+kappa*dir_k-mu_sig*ones(nn,1)))^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line Search %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     len=fminsearch(P,1)
len=1;
    sig_new=sig_new+len*dir_k;
    norm_grad=norm(J_total'*b_total);
    if iter==1
        norm_grad_init=norm_grad;
    end
    disp(['norm grad = ' num2str(norm_grad)]);
    iter=iter+1
    trisurf(tris,x,y,sig_new)
view(2)
axis off
shading interp
caxis([min(sigma_true) max(sigma_true)])
drawnow
end


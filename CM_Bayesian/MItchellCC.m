clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Continuum model EIT inversion %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('distmesh')
rng(123)            

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESHING
h = .075;         % ***

% PRIORS
npriors = 5;                        % # priors to generate
% bb = 1*1e1;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))
% aa = 1*1e0;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))
% cc = 4*1e1;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))

% ELECTRODES & PATTERNS
n_runs = 16;                  % # runs i.e. current patterns
n_meas_pts = 15;              % # measurement points
gauss_w = .5;                 % gaussian width(?) for current patterns

% STOCHASTIC
% n_bae_samps = 50;             % (unused) # bayesian samples (?) for ***

% RECONSTRUCTION
max_iter = 20;                % max # iterations for reconstruction

% PHYSICAL
nl = 1.;          % noise level [%]
beta = 1;         % ***
sig_0 = 1;        % sigma_0: (background?) conductivity
mu_sig = 0;       % mu_sigma: (mean?) conductivity

% MODEL
% fwd_shape = 'circle'; % 'circle', 'torso', or 'torso_subject'
fwd_shape = 'torso_subject';
% fwd_shape = 'torso_subject';
inv_shape = 'circle';

% GEOMETRY-BASED
if strcmp(fwd_shape,'torso')
    % Objects
    obj_c = {[-.40, .00], [.40, .00]};      % object centres
    obj_w = {.35, .35};                     % object widths
    obj_m = {-1.0, -1.0};                   % object magnitudes
    % Priors
    prior_corr_length = .5;       % (correlation?) length for priors
elseif strcmp(fwd_shape,'torso_subject')
    % Objects
    obj_c = {[-72, .00], [72, .00]};        % object centres
    obj_w = {63, 63};                       % object widths
    obj_m = {-1.0, -1.0};                   % object magnitudes
    % Priors
    prior_corr_length = 90;       % (correlation?) length for priors
else  % circle
    % Objects
    obj_c = {[-.45, -.35], [ .50, .20]};    % object centres
    obj_w = {.35, .35};                     % object widths
    obj_m = {1.0, -1.0};                    % object magnitudes
    % Priors
    prior_corr_length = .5;       % (correlation?) length for priors
end


%%%%%%%%%%%%%%%%%%%%%%%% Geometry & Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOMETRY & MESH

if strcmp(fwd_shape,'torso')
    [Nodes,Tris,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
    Nodes = Nodes./max(max(abs(Nodes)));    % temp: scaling to r_max=1
    tgeom = triangulation(Tris,Nodes);
elseif strcmp(fwd_shape,'torso_subject')
    [Nodes,Tris,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
    tgeom = triangulation(Tris,Nodes);
    mesh_file = "export_torso_param_mesh.txt";
    subject_file = "export_torso_param_test.txt";
    tgeom = generate_subject_mesh(tgeom,full_bdy,mesh_file,subject_file);
    Nodes = tgeom.Points;
else  % circle
    mesh_shape = @(p)dcircle(p,0,0,1);                                           
    pfix = [];
    figure(101)
    [Nodes,~] = distmesh2d(mesh_shape,@huniform,h,[-1,-1;1,1],pfix); % plots
    full_bdy = Nodes(:,1).^2+Nodes(:,2).^2>0.999;                   % boundary nodes
    tgeom = delaunayTriangulation(Nodes);
    Tris = tgeom.ConnectivityList;
end

% All nodes
x = Nodes(:,1); y = Nodes(:,2);                   % cartesian coords of nodes
theta = acos(x);                                  % polar angles of nodes
nn = length(Nodes);                               % # nodes

% All boundary nodes
bdy_indx = find(full_bdy);                        % indices of boundary nodes
nn_b = sum(full_bdy);                             % # nodes on boundary
x_bdy = x(full_bdy); y_bdy = y(full_bdy);
theta_bdy = atan2(y_bdy,x_bdy);       % polar angles of boundary nodes
[~,ss_b] = sort(theta_bdy);                       % indices when theta-sorted
theta_bdy_p = [theta_bdy-2*pi; ...
               theta_bdy; ...
               theta_bdy+2*pi];                   % periodic theta neumann

% Robin boundary condition
if strcmp(fwd_shape,'torso'),               xtol_rob = .1; 
elseif strcmp(fwd_shape,'torso_subject'),   xtol_rob = 18;     
else,                                       xtol_rob = .2;
end
bdy_rob = full_bdy & abs(x)<=xtol_rob & y<0;            % selection of robin BC
bdy_rob_indx = find(bdy_rob);                     % indices of robin BC
nn_rob = sum(bdy_rob);                            % # nodes with robin BC
[~,ss_rob] = sort(atan2(y(bdy_rob),x(bdy_rob)));  % indices when theta-sorted

% Neumann boundary condition
bdy_neu = full_bdy & ~bdy_rob;                    % neumann BC if not robin
bdy_neu_indx = find(bdy_neu);                     % indices of neumann BC
nn_neu = sum(bdy_neu);                            % # nodes with neumann BC
[~,ss_neu] = sort(atan2(y(bdy_neu),x(bdy_neu)));  % indices when theta-sorted

% ELECTRODES
theta_meas = linspace(0,2*pi,n_meas_pts+1)';         % polar theta array
theta_meas = theta_meas(1:end-1);                    % for meas points
if strcmp(fwd_shape,'torso') || strcmp(fwd_shape,'torso_subject')
    % interpolate between boundary nodes
    meas_pts = .99*[interp1(theta_bdy_p,[x_bdy; x_bdy; x_bdy],theta_meas-pi), ...
                    interp1(theta_bdy_p,[y_bdy; y_bdy; y_bdy],theta_meas-pi)];
else
    meas_pts = [.99*cos(theta_meas),.99*sin(theta_meas)];    % (x,y) coords
end


%%%%%%%%%%%%%%%%%%%%%%%%% Meas & Bdy Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neumann boundary condition   
neu_op = zeros(nn_neu,nn);                      % operator for neumann nodes
neu_op(1:nn_neu,bdy_neu_indx) = eye(nn_neu);    % 1s identify neumann nodes
theta_neu_c = atan2(y(bdy_neu),x(bdy_neu));     % theta of neumann nodes    (_c=circle?)
theta_neu_p = [theta_neu_c-2*pi; ...
               theta_neu_c; ...
               theta_neu_c+2*pi];               % periodic theta neumann

% Robin boundary condition
rob_op = zeros(nn_rob,nn);                      % operator for robin nodes
rob_op(1:nn_rob,bdy_rob_indx) = eye(nn_rob);    % 1s identify robin nodes

% Measurement electrodes
meas_op = zeros(n_meas_pts,nn);               % operator for meas nodes
for ii = 1:n_meas_pts 
    % find the triangle enclosing the measurement point    
    [ti,barycoords] = pointLocation(tgeom,meas_pts(ii,:));
    % store the nodal weights (barycentric) in the operator
    meas_op(ii,Tris(ti,:)) = barycoords;
end

%%%%%%%%%%% M_neu %%%%%%%%%%%
% (matrix?) for Neumann BC
[pM_n,~,i_Mn] = PreMassCurved1D(Nodes(bdy_neu_indx(ss_neu),:),0,0);         % ?
M_neu = sparse(i_Mn(:,1),i_Mn(:,2),pM_n*ones(nn_neu,1),nn_neu,nn_neu);
M_N = sparse([],[],0,nn,nn);
M_N(bdy_neu_indx(ss_neu),bdy_neu_indx(ss_neu)) = M_neu;

%%%%%%%%%%% M_rob %%%%%%%%%%%
% (matrix?) for Robin BC
[pM_r,~,i_Mr] = PreMassCurved1D(Nodes(bdy_rob_indx(ss_rob),:),0,0);             % ?
M_rob = sparse(i_Mr(:,1),i_Mr(:,2),pM_r*ones(nn_rob,1),nn_rob,nn_rob);
M_R = sparse([],[],0,nn,nn);
M_R(bdy_rob_indx(ss_rob),bdy_rob_indx(ss_rob)) = M_rob;

%%%%%%%%%%% M and K %%%%%%%%%%%
% (stiffness matrix) and (matrix?) for non-BC nodes
[pM,~,i_M] = Premass2D(Nodes,Tris,0);
M = sparse(i_M(:,1),i_M(:,2),pM*ones(nn,1),nn,nn);
[pK,SpK,i_K] = Prestiff2D(Nodes,Tris,1);
K = sparse(i_K(:,1),i_K(:,2),pK*ones(nn,1),nn,nn);

%%%%%%%%%%% F %%%%%%%%%%%
% Forcing matrix: generate current patterns for all neumann BCs
g = zeros(nn_neu,n_runs);
% (MP updated for domain continuity)
for ii = 1:n_runs
    % 1D gaussians along the boundary - g(theta_neu)
    % find the forcing function centres on the domain [-pi, pi]
    centre_pos = -pi + mod(  -pi+((ii-1)/n_runs)*2*pi, 2*pi);
    centre_neg = -pi + mod(-2*pi+((ii-1)/n_runs)*2*pi, 2*pi);
    % generate + and - currents on periodic domain [-3pi, 3pi]
    % function: phi=gaussian(centre,width,height,x)
    g_pos = gaussian(centre_pos,gauss_w,1,theta_neu_p);
    g_neg = gaussian(centre_neg,gauss_w,1,theta_neu_p);
    % sum the currents and collapse back to domain [-pi, pi]
    g(:,ii) = sum(reshape(g_pos-g_neg,[nn_neu,3]),2);
    % (optional) plot
    % figure(ii), plot(theta_neu_c(ss_neu), g(ss_neu,ii))
    % drawnow
end
clear centre_pos centre_neg g_pos g_neg
%  (original)
% for ii = 1:n_runs
%     g(:,ii) = gaussian(-pi+mod(-pi+(ii-1)/n_runs*2*pi,2*pi),gauss_w,1,theta_neu)-...
%         gaussian(-pi+mod(-2*pi+(ii-1)/n_runs*2*pi,2*pi),gauss_w,1,theta_neu);
% end
% Forcing matrix: (mass matrix?) operating on g
F = M_N*(neu_op'*(g));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian conductivity distribution (?)
G_sig = exp(-1/prior_corr_length^2*((x-x').^2+(y-y').^2)) + 1e-6*eye(nn);
% Lower(?) triangular of (?) sigma (Cholesky factorisation)
L_sig = chol(inv(G_sig));

% Plot the priors
rd = randn(nn,npriors);             % random for each node and prior
draw = mu_sig+L_sig\rd;             % to visualise the priors
subpltdim = ceil(sqrt(npriors));    % dimension of subplot to show priors
for ii = 1:npriors
    figure(102); subplot(subpltdim,subpltdim,ii)
    trisurf(Tris,x,y,draw(:,ii))
    view(2), shading interp, axis off, daspect([1,1,1])
    cb = colorbar; cb.TickLabelInterpreter = 'latex';
    caxis([min(draw(:)),max(draw(:))]), set(gca,'fontsize',14)
    title(strcat('prior sample',{' '},int2str(ii)),'interpreter','latex')
    drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Data synthesis s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Ground truth %%%%%%%%%%%
% true conductivity distribution
sig_t = mu_sig*ones(nn,1)+gaussian(obj_c{1},obj_w{1},obj_m{1},Nodes)+...
    gaussian(obj_c{2},obj_w{2},obj_m{2},Nodes);
% ?
K_t = sparse(i_K(:,1),i_K(:,2),pK*(sig_0*exp(sig_t)),nn,nn);
% ? --- WHY ROBIN MASS?
u = (K_t+beta*M_R)\F;

%%%%%%%%%%% Simulated data %%%%%%%%%%%
data_t = meas_op*u;                           % voltages at electrodes
size_data = size(data_t);                     % # electrode nodes
delta_noise = nl/100*(max(data_t(:)));        % noise amplitude
data = data_t+delta_noise*randn(size_data);   % noisy voltages
Le = 1/delta_noise;                           % ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Initialise reconstruction ? %%%%%%%%%%%
% sig_init = mu_sig*ones(nn,1);                 % initial guess of sigma (0's)
% sig_new = sig_init;                         % (unused)
data = data(:);

%%%%%%%%%%% (inverse geometry?) %%%%%%%%%%%
% Redefine variables for the assumed geometry (?)
% close all
% epsi = zeros(n_meas_pts*n_runs,n_bae_samps);  % (unused)
% x = x_circ;                                   % (undefined)
% y = y_circ;                                   % (undefined)
% Nodes = [x,y];                                % (redundant)
% theta_neu = atan2(y(bdy_neu),x(bdy_neu));     % (redundant)

% %%%%%%%%%%% Priors %%%%%%%%%%%
% % 
% [pM_prior,~,i_M_prior] = Premass2D(Nodes,Tris,0);
% M_prior = sparse(i_M_prior(:,1),i_M_prior(:,2),pM_prior*ones(nn,1),nn,nn);
% [pK_prior,~,i_K_prior] = Prestiff2D(Nodes,Tris,0);
% K_prior = sparse(i_K_prior(:,1),i_K_prior(:,2),pK_prior*ones(nn,1),nn,nn);
% L_prior = cc*(aa*K_prior+bb*M_prior);

%%%%%%%%%%% M_neu %%%%%%%%%%%
%
[pM_n,~,i_Mn] = PreMassCurved1D(Nodes(bdy_neu_indx(ss_neu),:),0,0);
M_neu = sparse(i_Mn(:,1),i_Mn(:,2),pM_n*ones(nn_neu,1),nn_neu,nn_neu);
M_N = sparse([],[],0,nn,nn);
M_N(bdy_neu_indx(ss_neu),bdy_neu_indx(ss_neu)) = M_neu;

%%%%%%%%%%% M_rob %%%%%%%%%%%
%
[pM_r,~,i_Mr] = PreMassCurved1D(Nodes(bdy_rob_indx(ss_rob),:),0,0);
M_rob = sparse(i_Mr(:,1),i_Mr(:,2),pM_r*ones(nn_rob,1),nn_rob,nn_rob);
M_R = sparse([],[],0,nn,nn);
M_R(bdy_rob_indx(ss_rob),bdy_rob_indx(ss_rob)) = M_rob;

%%%%%%%%%%% M, K_k, and S %%%%%%%%%%%
% repeat of above?
[pM,~,i_M] = Premass2D(Nodes,Tris,0);
M = sparse(i_M(:,1),i_M(:,2),pM*ones(nn,1),nn,nn);
[pK,~,i_K] = Prestiff2D(Nodes,Tris,0);
% below differs from original (K_k, and S)
K_k = sparse(i_K(:,1),i_K(:,2),pK*(sig_0*ones(nn,1)),nn,nn);
u_new = (K_k+beta*M_R)\F;
d_L = (meas_op/(K_k+beta*M_R));
S = sparse(n_runs*n_meas_pts,nn);
for jj = 1:n_runs
    S(1+(jj-1)*n_meas_pts:jj*n_meas_pts,:) = -d_L*reshape(SpK*u_new(:,jj),nn,nn);
end

%%%%%%%%%%% F %%%%%%%%%%%
%
g = zeros(nn_neu,n_runs);
% (MP updated for domain continuity)
for jj = 1:n_runs
    % 1D gaussians along the boundary - g(theta_neu)
    % find the forcing function centres on the domain [-pi, pi]
    centre_pos = -pi + mod(  -pi + ((jj-1)/n_runs)*2*pi, 2*pi);
    centre_neg = -pi + mod(-2*pi + ((jj-1)/n_runs)*2*pi, 2*pi);
    % generate + and - currents on periodic domain [-3pi, 3pi]
    % function: phi=gaussian(centre,width,height,x)
    g_pos = gaussian(centre_pos,gauss_w,1,theta_neu_p);
    g_neg = gaussian(centre_neg,gauss_w,1,theta_neu_p);
    % sum the currents and collapse back to domain [-pi, pi]
    g(:,jj) = sum(reshape(g_pos-g_neg,[nn_neu,3]),2);
    % (optional) plot
    % figure(ii), plot(theta_neu_c(ss_neu), g(ss_neu,ii))
    % drawnow
end
clear centre_pos centre_neg g_pos g_neg
% (original)
% for jj = 1:n_runs
%     g(:,jj) = gaussian(-pi+mod(-pi+(jj-1)/n_runs*2*pi,2*pi),gauss_w,1,theta_neu)-...
%         gaussian(-pi+mod(-2*pi+(jj-1)/n_runs*2*pi,2*pi),gauss_w,1,theta_neu);
% end
F = M_N*(neu_op'*(g));

%%
sig_new = mu_sig*ones(nn,1);
iter = 0;
norm_grad = 1;
norm_grad_init = 1;

% Plot the forward and inverse models
figure(104)
triplot(tgeom, 'Color', [.5 .5 .5]), hold on
axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
view(2), set(gca,'fontsize',14)
title('Forward model','interpreter','latex')
plot(x_bdy(ss_b),y_bdy(ss_b),'-k.')
plot(x(bdy_rob_indx),y(bdy_rob_indx),'or')
plot(x(bdy_neu_indx),y(bdy_neu_indx),'og')
plot(meas_pts(:,1),meas_pts(:,2),'xb','MarkerSize',12,'LineWidth',3)
hold off

% figure(103), subplot(2,2,1), triplot(tgeom, 'Color', [.5 .5 .5])
% axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
% view(2), set(gca,'fontsize',14)
% title('Forward model','interpreter','latex')

% Plot the ground truth
figure(103), subplot(2,2,3), trisurf(Tris,x,y,sig_t)
axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
view(2), shading interp, set(gca,'fontsize',14)
title('Ground truth','interpreter','latex')
% cb = colorbar; cb.TickLabelInterpreter = 'latex';
drawnow

% Plot the initial reconstruction
figure(103), subplot(2,2,4), trisurf(Tris,x,y,sig_new)
axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
view(2), shading interp, set(gca,'fontsize',14)
drawnow, pause(1)

% Iterate solution
while iter <= max_iter && norm_grad/norm_grad_init>1e-7
    
    K_k = sparse(i_K(:,1), i_K(:,2), pK*exp(sig_new), nn, nn);% (nn x nn) sparse
    u_new = (K_k+beta*M_R)\F;
    d_L = (meas_op/(K_k+beta*M_R));
    S = sparse(n_runs*n_meas_pts,nn);
    for jj = 1:n_runs
        S(1+(jj-1)*n_meas_pts:jj*n_meas_pts,:) = -d_L*reshape(SpK*u_new(:,jj),nn,nn);
    end
    J_total = [Le*S*diag(exp(sig_new));L_sig];                              % prior here
    b_total = [Le*(reshape(meas_op*u_new,n_runs*n_meas_pts,1)-data);...
        L_sig*(sig_new-mu_sig*ones(nn,1))];                                 % prior here
    dir_k = -J_total\b_total;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Cost Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = @(kappa) norm(Le*(reshape(meas_op*((sparse(i_K(:,1),i_K(:,2),pK*...
        exp(sig_new+kappa*dir_k))+beta*M_R)\F),n_runs*n_meas_pts,1)-data))^2+...
        norm(L_sig*(sig_new+kappa*dir_k-mu_sig*ones(nn,1)))^2;              % prior here

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line Search %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = fminsearch(P,1);
    sig_new = sig_new+len*dir_k;
    norm_grad = norm(J_total'*b_total);

    % Update iterations
    if iter == 1, norm_grad_init = norm_grad; end
    iter = iter+1;
    disp(['iter = ' num2str(iter) ',' char(9) ' norm grad = ' num2str(norm_grad)]);
    
    % Plot latest reconstruction
    figure(103), subplot(2,2,4), trisurf(Tris,x,y,sig_new)
    axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
    view(2), shading interp, set(gca,'fontsize',14)
    title('Latest iteration','interpreter','latex')
    drawnow, pause(1)
end







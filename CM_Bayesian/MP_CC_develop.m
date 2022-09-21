clear
clc
close all

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Continuum model EIT inversion %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('distmesh')
addpath('filexchange')
rng(123)            

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % MESHING
h = .075*150;         % ***

% PRIORS
npriors = 9;                        % # priors to generate
% bb = 1*1e1;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))
% aa = 1*1e0;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))
% cc = 4*1e1;       % coefficient (L_prior = cc*(aa*K_prior+bb*M_prior))

% ELECTRODES & PATTERNS
n_runs = 16;                  % # runs i.e. current patterns
n_meas_pts = 16;              % # measurement points
gauss_w = .2;%.5                 % gaussian width(?) for current patterns

% STOCHASTIC
% n_bae_samps = 50;             % (unused) # bayesian samples (?) for ***

% RECONSTRUCTION
max_iter = 20;                % max # iterations for reconstruction

% PHYSICAL
nl = 1.;          % noise level [%]
beta = 1;         % ***
sig_0 = 1;        % sigma_0: (?) conductivity
% mu_sig = 0;       % mu_sigma: (mean?) conductivity

% MODEL
% fwd_shape = 'circle'; % 'circle', 'torso', 'torso_subject',
% 'subject_truth' (in development)
fwd_shape = 'population_as_truth';
% fwd_shape = 'subject_truth';
% fwd_shape = 'torso_subject';
inv_shape = 'population';
% inv_shape = 'subject_predicted';

% % FORWARD GEOMETRY-BASED
% if strcmp(fwd_shape,'torso')
%     % Objects
%     obj_c = {[-.40, .00], [.40, .00]};      % object centres
%     obj_w = {.35, .35};                     % object widths
%     obj_m = {-1.0, -1.0};                   % object magnitudes
% elseif strcmp(fwd_shape,'torso_subject')
%     % Objects
%     obj_c = {[-72, .00], [72, .00]};        % object centres
%     obj_w = {63, 63};                       % object widths
%     obj_m = {-1.0, -1.0};                   % object magnitudes
% else  % circle
%     % Objects
%     obj_c = {[-.45, -.35], [ .50, .20]};    % object centres
%     obj_w = {.35, .35};                     % object widths
%     obj_m = {1.0, -1.0};                    % object magnitudes
% end

% INVERSE GEOMETRY-BASED
if strcmp(fwd_shape,'circle')
    prior_corr_length = .5;       % (correlation?) length for priors
else
    prior_corr_length = 1e2;       % (correlation?) length for priors
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packing parameters
param_m = {h};
param_e = {n_runs,n_meas_pts,gauss_w};
param_p = {nl,beta,sig_0}; %,mu_sig};
param_o = {0,0,0}; %{obj_c,obj_w,obj_m};
% param_m = {.075};
% param_e = {16,15,.5};
% param_p = {1.,1,1,0};
% param_o = {obj_c,obj_w,obj_m};
parameters = {param_m,param_e,param_p,param_o};

% test
% hiding the forward soltion to prevent any "cheating"
[data,sig_t_F,Le] = generate_forward_data(fwd_shape,parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Inverse Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOMETRY & MESH
[tgeom,full_bdy,var_r] = generate_mesh_inv(inv_shape);
Nodes = tgeom.Points;
Tris = tgeom.ConnectivityList;

% Define location of robin BC
% def_rob = {[1.45*pi 1.55*pi], };
def_rob = {[(0.45-1/8)*pi (0.55-1/8)*pi], };

% Generate the inverse model
geometry = {tgeom,Nodes,Tris,full_bdy};
[model_vars] = generate_model(geometry,{def_rob,n_meas_pts});

% Unpack all variables
[vars_all,vars_bdy,vars_rob,vars_neu,vars_mea,vars_op] = deal(model_vars{:});
[x,y,theta,nn] = deal(vars_all{:});
[x_bdy,y_bdy,theta_bdy,bdy_indx,nn_bdy,ss_b,~] = deal(vars_bdy{:});  % bdy_indx
[bdy_rob,bdy_rob_indx,nn_rob,ss_rob] = deal(vars_rob{:});
[bdy_neu,bdy_neu_indx,nn_neu,ss_neu,theta_neu_c,theta_neu_p] = deal(vars_neu{:});
[theta_meas,meas_pts] = deal(vars_mea{:});
[neu_op, rob_op, meas_op] = deal(vars_op{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create priors, mean shape info, and regularisation matrix (L_sig)
[mu_sig, L_sig] = generate_priors(nn,Tris,x,y,theta,npriors,prior_corr_length,true);
% stop


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
% Mass matrix for Neumann BC
% MP modified use of bdy_indicator
% test if the first and last boundary points are both in the neumann bc
% turn on boundary indicator if so
% ismember([bdy_indx(ss_b(1)) bdy_indx(ss_b(end))],bdy_neu_indx)
if all(ismember([bdy_indx(ss_b(1)) bdy_indx(ss_b(end))],bdy_neu_indx))
    [pM_n,~,i_Mn] = PreMassCurved1D(Nodes(bdy_neu_indx(ss_neu),:),1,0);
else
    [pM_n,~,i_Mn] = PreMassCurved1D(Nodes(bdy_neu_indx(ss_neu),:),0,0);
end
M_neu = sparse(i_Mn(:,1),i_Mn(:,2),pM_n*ones(nn_neu,1),nn_neu,nn_neu);
M_N = sparse([],[],0,nn,nn);
M_N(bdy_neu_indx(ss_neu),bdy_neu_indx(ss_neu)) = M_neu;

%%%%%%%%%%% M_rob %%%%%%%%%%%
% Mass matrix for Robin BC
% MP modified use of bdy_indicator
% test if the first and last boundary points are both in the neumann bc
% turn on boundary indicator if so
% ismember([bdy_indx(ss_b(1)) bdy_indx(ss_b(end))],bdy_rob_indx)
if all(ismember([bdy_indx(ss_b(1)) bdy_indx(ss_b(end))],bdy_rob_indx))
    [pM_r,~,i_Mr] = PreMassCurved1D(Nodes(bdy_rob_indx(ss_rob),:),1,0);
else
    [pM_r,~,i_Mr] = PreMassCurved1D(Nodes(bdy_rob_indx(ss_rob),:),0,0);
end
M_rob = sparse(i_Mr(:,1),i_Mr(:,2),pM_r*ones(nn_rob,1),nn_rob,nn_rob);
M_R = sparse([],[],0,nn,nn);
M_R(bdy_rob_indx(ss_rob),bdy_rob_indx(ss_rob)) = M_rob;

%%%%%%%%%%% M & K %%%%%%%%%%%
[pM,~,i_M] = Premass2D(Nodes,Tris,0);
M = sparse(i_M(:,1),i_M(:,2),pM*ones(nn,1),nn,nn);
[pK,SpK,i_K] = Prestiff2D(Nodes,Tris,1);  % [pK,~,i_K] = Prestiff2D(Nodes,Tris,0);
% below differs from original (K_k, and S)
K_k = sparse(i_K(:,1),i_K(:,2),pK*(sig_0.*ones(nn,1)),nn,nn);

%%%%%%%%%%% F %%%%%%%%%%%
%
g = zeros(nn_neu,n_runs);
for ii = 1:n_runs
    % 1D gaussians along the boundary - g(theta_neu)
    % find the forcing function centres on the domain [-pi, pi]
    % opposite electrode index
    oi = mod(ii + floorDiv(n_runs,2)-1,n_runs)+1;
    centre_pos = mod(theta_meas(ii)-pi,2*pi)-pi;
    centre_neg = mod(theta_meas(oi)-pi,2*pi)-pi;
    % generate + and - currents on periodic domain [-3pi, 3pi]
    % function: phi=gaussian(centre,width,height,x)
    g_pos = gaussian(centre_pos,gauss_w,1,theta_neu_p);
    g_neg = gaussian(centre_neg,gauss_w,1,theta_neu_p);
    % sum the currents and collapse back to domain [-pi, pi]
    g(:,ii) = sum(reshape(g_pos-g_neg,[nn_neu,3]),2);
    % % (optional) plot
    % figure(1e3+ii), plot(theta_neu_c(ss_neu), g(ss_neu,ii))
    % drawnow
end
clear centre_pos centre_neg g_pos g_neg
F = M_N*(neu_op'*(g));

% %%%%%%%%%%% S %%%%%%%%%%%
% u_new = (K_k+beta*M_R)\F;
% d_L = (meas_op/(K_k+beta*M_R));
% S = sparse(n_runs*n_meas_pts,nn);
% for jj = 1:n_runs
%     S(1+(jj-1)*n_meas_pts:jj*n_meas_pts,:) = -d_L*reshape(SpK*u_new(:,jj),nn,nn);
% end


%% Initiate solution
sig_new = mu_sig; %*ones(nn,1);
iter = 0;
norm_grad = 1;
norm_grad_init = 1;

% Update plots 
%  Plot the inverse model
% figure(104)
figure(103), subplot(2,2,2)
triplot(tgeom, 'Color', [.5 .5 .5]), hold on
axis off, daspect([1,1,1])
view(2), set(gca,'fontsize',10)
title('Inverse model','interpreter','latex')
plot(x_bdy(ss_b),y_bdy(ss_b),'-k.')
plot(x(bdy_rob_indx),y(bdy_rob_indx),'or')
plot(x(bdy_neu_indx),y(bdy_neu_indx),'og')
plot(meas_pts(:,1),meas_pts(:,2),'xb','MarkerSize',12,'LineWidth',3)
for i = 1:length(meas_pts)
    text(1.2*meas_pts(i,1),1.2*meas_pts(i,2),num2str(i),...
        'Color','b','FontSize',10)
end
hold off

% Plot the initial reconstruction
figure(103), subplot(2,2,4), trisurf(Tris,x,y,sig_new)
axis off, daspect([1,1,1]), %caxis([min(sig_t) max(sig_t)])
% xlim([-250 250]), ylim([-200 200])
view(2), shading interp, set(gca,'fontsize',10)
drawnow, pause(1)


%% Iterate solution
tol = 1e-4;  % 1e-7
while iter <= max_iter && norm_grad/norm_grad_init>tol
    
    K_k = sparse(i_K(:,1), i_K(:,2), pK*exp(sig_new), nn, nn);% (nn x nn) sparse
    u_new = (K_k+beta*M_R)\F;
    d_L = (meas_op/(K_k+beta*M_R));
    S = sparse(n_runs*n_meas_pts,nn);
    for jj = 1:n_runs
        S(1+(jj-1)*n_meas_pts:jj*n_meas_pts,:) = -d_L*reshape(SpK*u_new(:,jj),nn,nn);
    end
    J_total = [Le*S*diag(exp(sig_new));L_sig];                              % prior here
    b_total = [Le*(reshape(meas_op*u_new,n_runs*n_meas_pts,1)-data);...
        L_sig*(sig_new-mu_sig)]; %*ones(nn,1))];                                 % prior here
    dir_k = -J_total\b_total;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Cost Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = @(kappa) norm(Le*(reshape(meas_op*((sparse(i_K(:,1),i_K(:,2),pK*...
        exp(sig_new+kappa*dir_k))+beta*M_R)\F),n_runs*n_meas_pts,1)-data))^2+...
        norm(L_sig*(sig_new+kappa*dir_k-mu_sig))^2; %*ones(nn,1)))^2;              % prior here

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
    axis off, daspect([1,1,1]), %caxis([min(sig_t) max(sig_t)])
%     xlim([-250 250]), ylim([-200 200])
    view(2), shading interp, set(gca,'fontsize',10)
    title('Latest iteration','interpreter','latex')
    drawnow, pause(1)
end
figure(103), subplot(2,2,4), title('Final iteration','interpreter','latex')


%% Posterior
%  Posterior calculations
G_post=inv(J_total'*J_total);               % 
SD_post=diag(G_post).^(1/2);                % standard deviations
L_post=chol(G_post,'lower');
n_post_samples = 100;
sig_samples=sig_new+L_post*randn(nn,n_post_samples);

%  Mahalanobis distance of solution
sig_t_invmesh = sig_t_F(x,y);
m_dist = (((sig_t_invmesh-sig_new)')/(G_post)*(sig_t_invmesh-sig_new)).^0.5;
fprintf("\nPrior correlation length: d = %9.2e\n", prior_corr_length)
fprintf("\nMahalanobis distance:\t  d = %7.4f\n", m_dist)
% d^2 follows the chi-squared distribution with n degrees of freedom, 
% where n is the number of dimensions of the normal distribution.
%  n =745?

%  Total variation (?)
%  first calculate mesh element areas
Tri_areas = zeros([length(Tris) 1]);
for i = 1:length(Tri_areas)
    x1 = Nodes(Tris(i,1),1); y1 = Nodes(Tris(i,1),2);
    x2 = Nodes(Tris(i,2),1); y2 = Nodes(Tris(i,2),2);
    x3 = Nodes(Tris(i,3),1); y3 = Nodes(Tris(i,3),2);
    Tri_areas(i) = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
end

%  then gradient and total variance at each element (relative to mean) 
%  for each posterior sample
tot_vars = zeros([length(n_post_samples) 1]);
for ii=1:n_post_samples
    [dFx,dFy] = trigradient(Tris,x,y,sig_samples(:,ii)-sig_new,'face'); % at faces!
    tot_vars(ii) = sum((abs(dFx)+abs(dFy)).*Tri_areas)/sum(Tri_areas);
end
fprintf("\nTotal variation:     \t  d = %7.4f \x00B1 %7.4f\n\n", mean(tot_vars), std(tot_vars))


%% Posterior plotting
scrnwh = get(0,'screensize');
f5w = 900; f5h = 400;
cax_lims1 = [min(sig_t_invmesh) max(sig_t_invmesh)];
cax_lims2 = [min(min(sig_samples)) max(max(sig_samples))];

% % prior mean
figure(104); subplot(3,2,3), trisurf(Tris,x,y,mu_sig) %truth needs own tris,x,y...
axis equal, axis off, caxis(cax_lims1)
% xlim([-250 250]), ylim([-200 200])
view(2), shading interp
title('Prior Mean','interpreter','latex')

% % reconstructed mean
% sig_t_LS = sig_t; %get_level_set(sig_t);
fig = figure(104); subplot(3,2,2), trisurf(Tris,x,y,sig_new) 
axis equal, axis off, caxis(cax_lims1)
title('Reconstructed Mean','interpreter','latex')
% xlim([-250 250]), ylim([-200 200])
view(2), shading interp

% % reconstructed stdev
sig_t_LS = sig_t_invmesh; %get_level_set(sig_t_invmesh);
figure(104); subplot(3,2,4), trisurf(Tris,x,y,SD_post) 
axis equal, axis off, caxis(cax_lims1)
title('Posterior Std. Dev.','interpreter','latex')
% xlim([-250 250]), ylim([-200 200])
view(2), shading interp

% % solution slice
x_samp = linspace(min(x),max(x),1e3);
y_samp1 = -30*ones(size(x_samp));
y_samp2 =  30*ones(size(x_samp));
F = scatteredInterpolant(x,y,sig_new); z_samp1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,SD_post); sd_samp1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,sig_samples(:,1)); zs1_samp1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,sig_samples(:,2)); zs2_samp1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,sig_samples(:,3)); zs3_samp1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,sig_t_invmesh); z_true1 = F(x_samp,y_samp1);
F = scatteredInterpolant(x,y,sig_new); z_samp2 = F(x_samp,y_samp2);
F = scatteredInterpolant(x,y,SD_post); sd_samp2 = F(x_samp,y_samp2);
F = scatteredInterpolant(x,y,sig_samples(:,1)); zs1_samp2 = F(x_samp,y_samp2);
F = scatteredInterpolant(x,y,sig_samples(:,2)); zs2_samp2 = F(x_samp,y_samp2);
F = scatteredInterpolant(x,y,sig_samples(:,3)); zs3_samp2 = F(x_samp,y_samp2);
F = scatteredInterpolant(x,y,sig_t_invmesh); z_true2 = F(x_samp,y_samp2);
z_samp1_p99 = z_samp1+2.576*sd_samp1;
z_samp1_m99 = z_samp1-2.576*sd_samp1;
z_samp2_p99 = z_samp2+2.576*sd_samp2;
z_samp2_m99 = z_samp2-2.576*sd_samp2;
figure(105); subplot(3,2,1:2), trisurf(Tris,x,y,sig_new), hold on
plot3(x_samp,y_samp1,2*abs(z_samp1),'k-')
plot3(x_samp,y_samp2,2*abs(z_samp2),'k-')
text(max(x_samp),max(y_samp1),max(2*abs(z_samp1)),'A',...
     'HorizontalAlignment','left','VerticalAlignment','middle')
text(max(x_samp),max(y_samp2),max(2*abs(z_samp2)),'B',...
     'HorizontalAlignment','left','VerticalAlignment','middle')
set(gcf,'Position',[scrnwh(3)-f5w scrnwh(4)-f5h f5w f5h])
axis equal, axis off, caxis(cax_lims1)
title('MAP Estimate','interpreter','latex')
view(2), shading interp, hold off
figure(105); subplot(3,2,[3 5]), hold on
fill([x_samp fliplr(x_samp)], [z_samp1_p99, fliplr(z_samp1_m99)], ...
     .8*[1 1 1], 'LineStyle','none', 'FaceAlpha',0.5);
plot(x_samp,zs1_samp1,':','LineWidth',2,'Color',.5*[1 1 1]);
plot(x_samp,zs2_samp1,':','LineWidth',2,'Color',.5*[1 1 1]);
plot(x_samp,zs3_samp1,':','LineWidth',2,'Color',.5*[1 1 1]);
% plot(x_samp,z_samp1_p99,'--','LineWidth',2,'Color',1*[1 1 1]);
% plot(x_samp,z_samp1_m99,'--','LineWidth',2,'Color',1*[1 1 1]);
% area(x_samp,z_true1,'FaceColor', [.7 .7 .7]);
plot(x_samp,z_true1,'r-','LineWidth',2);
plot(x_samp,z_samp1,'-','LineWidth',2,'Color',.0*[1 1 1]);
xlim([min(x_samp) max(x_samp)]), xlabel('x'), ylabel('${\sigma}$')
title('Section A','interpreter','latex')
plots=get(gca, 'Children');
legend([plots(2) plots(1) plots(6) plots(5)], ...
       {'True','MAP','${99\%}$ CI','Samples'}, ...
       'Location','SouthOutside','Orientation','horizontal')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'layer', 'top')
figure(105); subplot(3,2,[4 6]), hold on
fill([x_samp fliplr(x_samp)], [z_samp2_p99, fliplr(z_samp2_m99)], ...
     .8*[1 1 1], 'LineStyle','none', 'FaceAlpha',0.3);
plot(x_samp,zs1_samp2,':','LineWidth',2,'Color',.5*[1 1 1]);
plot(x_samp,zs2_samp2,':','LineWidth',2,'Color',.5*[1 1 1]);
plot(x_samp,zs3_samp2,':','LineWidth',2,'Color',.5*[1 1 1]);
% plot(x_samp,z_samp2_p99,'--','LineWidth',2,'Color',1*[1 1 1]);
% plot(x_samp,z_samp2_m99,'--','LineWidth',2,'Color',1*[1 1 1]);
% area(x_samp,z_true2,'FaceColor', [.7 .7 .7]);
plot(x_samp,z_true2,'r-','LineWidth',2);
plot(x_samp,z_samp2,'-','LineWidth',2,'Color',.0*[1 1 1]);
xlim([min(x_samp) max(x_samp)]), xlabel('x'), ylabel('${\sigma}$')
title('Section B','interpreter','latex')
plots=get(gca, 'Children');
legend([plots(2) plots(1) plots(6) plots(5)], ...
       {'True','MAP','${99\%}$ CI','Samples'}, ...
       'Location','SouthOutside','Orientation','horizontal')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'layer', 'top')
hold off


% % % posterior samples
% for ii=1:n_post_samples
%     sig_sample_LS = sig_samples(:,ii);%get_level_set(sig_samples(:,ii));
%     figure(104), subplot(3,2,6), trisurf(Tris,x,y,sig_sample_LS)
%     axis equal, axis off, caxis(cax_lims1)
%     title('Posterior Samples','interpreter','latex')
%     view(2), shading interp
%     drawnow
% 
%     % Capture the plot as an image 
%     frame = getframe(fig); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if ii == 1, imwrite(imind,cm,'eit_animated.gif','gif', 'Loopcount',inf); 
%     else, imwrite(imind,cm,'eit_animated.gif','gif','WriteMode','append'); end 
% 
% end

%% FUNCTIONS

function [data_LS] = get_level_set(data)
    data_LS=0*data;
    data_LS(data>.5)=1;
    data_LS(data<-.5)=-1;
end

function [tgeom,full_bdy,var_r] = generate_mesh_inv(inv_shape)

    if strcmp(inv_shape,'torso')
        [tgeom,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
        var_r = [];

    elseif strcmp(inv_shape,'torso_subject')
        [tgeom,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
        mesh_file = "./geom/bspfile_torso_mean.txt";
        subject_file = "./geom/bspfile_torso_sbj_test.txt";
        [tgeom,var_r] = generate_subject_mesh(tgeom,full_bdy,mesh_file,subject_file);
    
    elseif strcmp(inv_shape,'subject_predicted')
        trimesh_file = './geom/TEST_population_sample_mean/trimesh_0010';
        [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);
        bspfile_mesh = "./geom/TEST_population_sample_mean/bspfile_torso.txt";
        bspfile_subject = "./geom/TEST_HLA-H11303_predicted/bspfile_torso.txt";
        [tgeom,var_r] = generate_subject_mesh(tgeom,full_bdy,bspfile_mesh,bspfile_subject); 
    
    elseif strcmp(inv_shape,'population')
        trimesh_file = './geom/TEST_population_sample_mean/trimesh_0010';
        [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);
        var_r = [];
    
    else  % circle
        mesh_shape = @(p)dcircle(p,0,0,150);                                           
        pfix = [];
        figure(101)
        [Nodes,~] = distmesh2d(mesh_shape,@huniform,h,[-1,-1;1,1]*150,pfix); % plots
        full_bdy = (Nodes(:,1).^2+Nodes(:,2).^2).^.5>0.999*150;                   % boundary nodes
        tgeom = delaunayTriangulation(Nodes);
    
    end
end


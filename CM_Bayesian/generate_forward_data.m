function [data,sig_t_F,Le] = generate_forward_data(fwd_shape,parameters)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Continuum model EIT inversion %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('distmesh')
addpath('bspline')
rng(123)            

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters = {param_m,param_e,param_p,param_o}
[param_m,param_e,param_p,param_o] = deal(parameters{:});
h = param_m{1};
[n_runs,n_meas_pts,gauss_w] = deal(param_e{:});
[nl,beta,sig_0] = deal(param_p{:});
[obj_c,obj_w,obj_m] = deal(param_o{:});

% Define location of robin BC
% def_rob = {[1.45*pi 1.55*pi], };
def_rob = {[(0.45-1/8)*pi (0.55-1/8)*pi], };

%%%%%%%%%%%%%%%%%%%%%%%% Geometry & Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOMETRY & MESH
[tgeom, full_bdy] = generate_mesh_fwd(fwd_shape);
Nodes = tgeom.Points;
Tris = tgeom.ConnectivityList;  

% Generate the forward model
geometry = {tgeom,Nodes,Tris,full_bdy};
[model_vars] = generate_model(geometry,{def_rob,n_meas_pts});

% Unpack all forward model variables
[vars_all,vars_bdy,vars_rob,vars_neu,vars_mea,vars_op] = deal(model_vars{:});
[x,y,theta,nn] = deal(vars_all{:});
[x_bdy,y_bdy,theta_bdy,bdy_indx,~,ss_b,~] = deal(vars_bdy{:});  % bdy_indx
[~,bdy_rob_indx,nn_rob,ss_rob] = deal(vars_rob{:});
[~,bdy_neu_indx,nn_neu,ss_neu,theta_neu_c,theta_neu_p] = deal(vars_neu{:});
[theta_meas,meas_pts] = deal(vars_mea{:});
[neu_op,~,meas_op] = deal(vars_op{:});


%%%%%%%%%%%%%%%%%%%%%%%%% Meas & Bdy Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%% M and K %%%%%%%%%%%
% Mass and stiffness matrices for non-BC nodes
[pM,~,i_M] = Premass2D(Nodes,Tris,0);
M = sparse(i_M(:,1),i_M(:,2),pM*ones(nn,1),nn,nn);
[pK,~,i_K] = Prestiff2D(Nodes,Tris,0);  % [pK,SpK,i_K] = Prestiff2D(Nodes,Tris,1);
K = sparse(i_K(:,1),i_K(:,2),pK*ones(nn,1),nn,nn);

%%%%%%%%%%% F %%%%%%%%%%%
% Forcing matrix: generate current patterns for all neumann BCs
g = zeros(nn_neu,n_runs);
% (MP updated for domain continuity)
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
    % (optional) plot
    % figure(1e3+ii), plot(theta_neu_c(ss_neu), g(ss_neu,ii))
    % drawnow
end
clear centre_pos centre_neg g_pos g_neg
% Forcing matrix: (mass matrix?) operating on g
F = M_N*(neu_op'*(g));


% figure(999)
% theta_neu = theta(bdy_neu_indx); theta_rob = theta(bdy_rob_indx);
% for i = 1:16
% %     x_plot = x(full_bdy); y_plot = y(full_bdy); F_plot = F(full_bdy);
% %     subplot(4,4,i), plot3(x_plot(ss_b),y_plot(ss_b),F_plot(ss_b))
% %     axis off, daspect([1,1,1]), caxis([min(F(:)) max(F(:))])
% %     view(2), shading interp
%     F_bdy = F(full_bdy,i); F_neu = F(bdy_neu_indx,i); F_rob = F(bdy_rob_indx,i);
%     subplot(4,4,i), plot(theta_bdy(ss_b),F_bdy(ss_b),'g'), hold on
%     plot(theta_neu(ss_neu),F_neu(ss_neu),'.b','MarkerSize',3)
%     plot(theta_rob(ss_rob),F_rob(ss_rob),'.r','MarkerSize',3)
%     set(gca,'fontsize',8)
%     ylim([min(F(:)) max(F(:))])
%     xlabel('$\theta$','interpreter','latex')
%     ylabel(strcat('F($\theta$,',string(i),')'),'interpreter','latex')
%     hold off
% end

% figure(997)
% for i = 1:16
%     subplot(4,4,i), plot(theta_neu(ss_neu),g(ss_neu,i),'g'), hold on
%     g_all = neu_op'*(g);
%     plot(theta,g_all(:,i),'.b','MarkerSize',3)
%     set(gca,'fontsize',8)
%     ylim([min(g(:)) max(g(:))])
%     xlabel('$\theta$','interpreter','latex')
%     ylabel(strcat('F($\theta$,',string(i),')'),'interpreter','latex')
%     hold off
% end

% figure(998)
% for i = 1:16
% %     x_plot = x(full_bdy); y_plot = y(full_bdy); F_plot = F(full_bdy);
% %     subplot(4,4,i), plot3(x_plot(ss_b),y_plot(ss_b),F_plot(ss_b))
% %     axis off, daspect([1,1,1]), caxis([min(F(:)) max(F(:))])
% %     view(2), shading interp
%     x_neu = x(bdy_neu_indx); y_neu = y(bdy_neu_indx); 
%     x_rob = x(bdy_rob_indx); y_rob = y(bdy_rob_indx);
%     F_plot = F(full_bdy,i);
%     subplot(4,4,i), plot3(x_bdy(ss_b),y_bdy(ss_b),F_plot(ss_b),'g'), hold on
%     plot3(x_neu(ss_neu),y_neu(ss_neu),F_plot(ss_neu),'.b','MarkerSize',10)
%     plot3(x_rob(ss_rob),y_rob(ss_rob),F_plot(ss_rob),'.r','MarkerSize',10)
%     view(2), set(gca,'fontsize',8)
%     %ylim([min(F(:)) max(F(:))])
% %     xlabel('$\theta$','interpreter','latex')
% %     ylabel(strcat('F($\theta$,',string(i),')'),'interpreter','latex')
%     hold off
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Data synthesis s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Ground truth %%%%%%%%%%%
% true conductivity distribution
% sig_t = mu_sig*ones(nn,1) + gaussian(obj_c{1},obj_w{1},obj_m{1},Nodes)+...
%                             gaussian(obj_c{2},obj_w{2},obj_m{2},Nodes);

% % from b-spline...
% [data_left,shift_left,~] = read_bspline_file('export_lung_left_param_test.txt');
% [data_rght,shift_rght,~] = read_bspline_file('export_lung_rght_param_test.txt');
% % rho for the theta of each boundary node
% theta_lung = linspace(0,2*pi,201)';
% theta_lung = theta_lung(1:end-1);
% rtheta_left = evaluate_mesh_bspline(data_left(:,2), data_left(:,1), theta_lung);
% rtheta_rght = evaluate_mesh_bspline(data_rght(:,2), data_rght(:,1), theta_lung);
% [x_left,y_left] = pol2cart(theta_lung,rtheta_left);
% [x_rght,y_rght] = pol2cart(theta_lung,rtheta_rght);
% x_left = -x_left - shift_left(1); y_left = -y_left - shift_left(2);
% x_rght = -x_rght - shift_rght(1); y_rght = -y_rght - shift_rght(2);
% x_left = x_left - 165; y_left = y_left - 180;
% x_rght = x_rght - 165; y_rght = y_rght - 180;
% lung_left = inpolygon(x,y,x_left,y_left);
% lung_rght = inpolygon(x,y,x_rght,y_rght);
% % lung_left(y>mean(y)) = 0;  % add diseased state to lung: LEFT,  POSTERIOR
% % lung_left(y<mean(y)) = 0;  % add diseased state to lung: LEFT, ANTERIOR
% % lung_rght(y>mean(y)) = 0;  % add diseased state to lung: RIGHT, POSTERIOR
% % lung_rght(y<mean(y)) = 0;  % add diseased state to lung: RIGHT, ANTERIOR
% sig_t = mu_sig*ones(nn,1) + lung_left + lung_rght;

% from mesh seeds
xy_left = readmatrix('./geom/TEST_population_sample_mean/meshseed_lung_left_0005.csv');
xy_rght = readmatrix('./geom/TEST_population_sample_mean/meshseed_lung_right_0005.csv');
% xy_left = readmatrix('./geom/TEST_HLA-H11303_truth/meshseed_lung_left_0005.csv');
% xy_rght = readmatrix('./geom/TEST_HLA-H11303_truth/meshseed_lung_right_0005.csv');
lung_left = inpolygon(x,y,xy_left(:,1),xy_left(:,2));
lung_rght = inpolygon(x,y,xy_rght(:,1),xy_rght(:,2));
% lung_left(y>mean(y)) = 0;  % add diseased state to lung: LEFT,  POSTERIOR
% lung_left(y<mean(y)) = 0;  % add diseased state to lung: LEFT, ANTERIOR
lung_rght(y>mean(y)) = 0;  % add diseased state to lung: RIGHT, POSTERIOR
% lung_rght(y<mean(y)) = 0;  % add diseased state to lung: RIGHT, ANTERIOR
sig_t = 0*ones(nn,1) + lung_left + lung_rght;
sig_t(sig_t>1) = 1;
sig_t_F = scatteredInterpolant(x,y,sig_t);

% figure(999)
% plot(theta_lung,rtheta_left), hold on
% plot(theta_lung,rtheta_rght)
% plot(x_left,y_left), hold on
% plot(x_rght,y_rght)

% ?
K_t = sparse(i_K(:,1),i_K(:,2),pK*(sig_0.*exp(sig_t)),nn,nn);
% ? 
u = (K_t+beta*M_R)\F;

%%%%%%%%%%% Simulated data %%%%%%%%%%%
data_t = meas_op*u;                           % voltages at electrodes
size_data = size(data_t);                     % # electrode nodes
delta_noise = nl/100*(max(data_t(:)));        % noise amplitude
data = data_t+delta_noise*randn(size_data);   % noisy voltages
Le = 1/delta_noise;                           % ?


%%%%%%%%%%% Plot %%%%%%%%%%%
scrnwh = get(0,'screensize');
f3w = 900; f3h = 600; 
f4w = 700; f4h = 790; 

% Plot the forward model
% figure(104)
figure(103), subplot(2,2,1)
set(gcf,'Position',[f4w scrnwh(4)-f3h f3w f3h])
triplot(tgeom, 'Color', [.5 .5 .5]), hold on
axis off, daspect([1,1,1])
view(2), set(gca,'fontsize',10)
title('Forward model','interpreter','latex')
plot(x_bdy(ss_b),y_bdy(ss_b),'-k.')
plot(x(bdy_rob_indx),y(bdy_rob_indx),'or')
plot(x(bdy_neu_indx),y(bdy_neu_indx),'og')
plot(meas_pts(:,1),meas_pts(:,2),'xb','MarkerSize',12,'LineWidth',3)
for i = 1:length(meas_pts)
    text(1.2*meas_pts(i,1),1.2*meas_pts(i,2),num2str(i),...
        'Color','b','FontSize',10)
end
plot(xy_left(:,1),xy_left(:,2),'k'), plot(xy_rght(:,1),xy_rght(:,2),'k')
hold off

% Plot the ground truth
figure(103), subplot(2,2,3), trisurf(Tris,x,y,sig_t)
axis off, daspect([1,1,1]), caxis([min(sig_t) max(sig_t)])
% xlim([-250 250]), ylim([-200 200])
view(2), shading interp, set(gca,'fontsize',10)
title('Ground truth','interpreter','latex')
% cb = colorbar; cb.TickLabelInterpreter = 'latex';
drawnow

figure(104); subplot(3,2,1), trisurf(Tris,x,y,sig_t) %truth needs own tris,x,y...
axis equal, axis off, caxis([min(sig_t) max(sig_t)])
title('Ground Truth','interpreter','latex')
set(gcf,'Position',[0 scrnwh(4)-f4h f4w f4h])

% xlim([-250 250]), ylim([-200 200])
view(2), shading interp

end


function [tgeom,full_bdy] = generate_mesh_fwd(fwd_shape)

    if strcmp(fwd_shape,'torso')
        [tgeom,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
%         [tgeom,full_bdy] = read_2D_tri_mesh('./geom/torso_mean_mesh_0010',0.7,false,false);
      
    elseif strcmp(fwd_shape,'torso_subject')
        [tgeom,full_bdy] = read_2D_stl_mesh('./geom/ST4_A',0.7,true,false);
        mesh_file = "export_torso_param_mesh.txt";
        subject_file = "export_torso_param_test.txt";
        [tgeom,~] = generate_subject_mesh(tgeom,full_bdy,mesh_file,subject_file);

    elseif strcmp(fwd_shape,'subject_truth')
        trimesh_file = './geom/TEST_HLA-H11303_truth/trimesh_0005';
        [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);

    elseif strcmp(fwd_shape,'prediction_as_truth')
        trimesh_file = './geom/TEST_population_sample_mean/trimesh_lung_0005';
        [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);
        bspfile_mesh = "./geom/TEST_population_sample_mean/bspfile_torso.txt";
        bspfile_subject = "./geom/TEST_HLA-H11303_predicted/bspfile_torso.txt";
        [tgeom,~] = generate_subject_mesh(tgeom,full_bdy,bspfile_mesh,bspfile_subject); 
    
    elseif strcmp(fwd_shape,'population_as_truth')
        trimesh_file = './geom/TEST_population_sample_mean/trimesh_lung_0005';
        [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);

    else  % circle
        mesh_shape = @(p)dcircle(p,0,0,150);                                           
        pfix = [];
        figure(101)
        [Nodes,~] = distmesh2d(mesh_shape,@huniform,h,[-1,-1;1,1],pfix); % plots
        full_bdy = Nodes(:,1).^2+Nodes(:,2).^2>0.999;                   % boundary nodes
        tgeom = delaunayTriangulation(Nodes);
    end
end
function [data,sig_t_F,Le] = generate_forward_data(parameters)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

addpath('distmesh','cem_bayesian')
rng(123)            

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param_m,param_e,param_d] = deal(parameters{:});
[fwd_shape] = deal(param_m{:});
[n_elec,n_patt,z_elec,w_elec,N_,MeasPattern] = deal(param_e{:});
[nl,sig_0] = deal(param_d{:});


%%%%%%%%%%%%%%%%%%%%%%%% Geometry & Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [tgeom, full_bdy] = generate_mesh_fwd(fwd_shape);
[tgeom,full_bdy] = read_2D_tri_mesh(fwd_shape,nan,false,false);
nodes = tgeom.Points;
tris = tgeom.ConnectivityList;  

%%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the forward model
params_m = {tgeom,nodes,tris,full_bdy};
params_e = {n_elec, w_elec};
fwd_params = {params_m,params_e};
[model_vars] = generate_model(fwd_params);

% Unpack all forward model variables
[vars_all,vars_bdy,vars_elec] = deal(model_vars{:});
[x,y,~,nn] = deal(vars_all{:});
[x_bdy,y_bdy,~,~,bdy_indx,ss_b,bdy_elems] = deal(vars_bdy{:});  % bdy_indx
[elec_pts, in_elec, n_per_elec] = deal(vars_elec{:});

% %TMEP
% plot_model('Forward model',123,0,tgeom,vars_bdy,vars_elec), hold on
% stop


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ground truth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Defined manually %%%%%%%%%%%%
%{obj_c,obj_w,obj_m};
% sig_t = mu_sig*ones(nn,1) + gaussian(obj_c{1},obj_w{1},obj_m{1},Nodes)+...
%                             gaussian(obj_c{2},obj_w{2},obj_m{2},Nodes);

%%%%%%%%%%% Defined from mesh seeds %%%%%%%%%%%
xy_left = -readmatrix('/hpc/mpag253/EIT/CM_Bayesian/geom/pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A/truth_H5977/mesh_seeds/mesh_seeds_lung_left_0005.csv');
xy_rght = -readmatrix('/hpc/mpag253/EIT/CM_Bayesian/geom/pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A/truth_H5977/mesh_seeds/mesh_seeds_lung_right_0005.csv');
% xy_left = readmatrix('./geom/TEST_HLA-H11303_truth/meshseed_lung_left_0005.csv');
% xy_rght = readmatrix('./geom/TEST_HLA-H11303_truth/meshseed_lung_right_0005.csv');
lung_left = inpolygon(x,y,xy_left(:,1),xy_left(:,2));
lung_rght = inpolygon(x,y,xy_rght(:,1),xy_rght(:,2));
% lung_left(y>mean(y)) = 0;  % add diseased state to lung: LEFT,  POSTERIOR
% lung_left(y<mean(y)) = 0;  % add diseased state to lung: LEFT, ANTERIOR
lung_rght(y>mean(y)) = 0;  % add diseased state to lung: RIGHT, POSTERIOR
% lung_rght(y<mean(y)) = 0;  % add diseased state to lung: RIGHT, ANTERIOR
sig_true = 0*ones(nn,1) + lung_left + lung_rght;
sig_true(sig_true>1) = 1;
sig_t_F = scatteredInterpolant(x,y,sig_true);  

% % Plot
% figure(999)
% plot(theta_lung,rtheta_left), hold on
% plot(theta_lung,rtheta_rght)
% plot(x_left,y_left), hold on
% plot(x_rght,y_rght)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate FEM model parameters
params_m = {nodes,tris,nn,bdy_indx,bdy_elems};
params_e = {n_elec,in_elec,n_per_elec,z_elec,N_};
params_p = {sig_true};
fem_params = {params_m,params_e,params_p};
[K,pK,SpK,i_K,M_B,C,D] = generate_fem_model(fem_params); 

% Solve
B = K+1/z_elec*M_B;
F = [zeros(nn,n_patt); N_'*MeasPattern];
A = [B,C;C',D];
uu = A\F;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Meas_rest = sparse((1:n_elec-1),(nn+1:nn+n_elec-1),1);
data_t = MeasPattern*N_*Meas_rest*uu;           % voltages at electrodes
size_data = size(data_t);                       % # electrode nodes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nl ~= 0
    delta_noise = nl/100*(max(data_t(:)));          % noise amplitude
    data = data_t+delta_noise*randn(size_data);     % noisy voltages
    Le = (1/delta_noise)*speye(n_patt*n_elec);      % (?)
else
    data = data_t;                                  % non-noisy voltages
    Le = 1e4*speye(n_patt*n_elec);                      % (?) VALIDATE!
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the forward model
plot_model('Forward model',123,0,tgeom,vars_bdy,vars_elec), hold on
plot([xy_left(:,1);xy_left(1,1)],[xy_left(:,2);xy_left(1,2)],'k')
plot([xy_rght(:,1);xy_rght(1,1)],[xy_rght(:,2);xy_rght(1,2)],'k'), hold off

% Plot the ground truth
plot_sigma('Ground truth',103,[1,2,1],tris,x,y,sig_true,0)
plot_sigma('Ground truth',104,[3,2,1],tris,x,y,sig_true,0)

end


% function [tgeom,full_bdy] = generate_mesh_fwd(fwd_shape)
% 
%     %trimesh_file = '/hpc/mpag253/EIT/CEM_Bayesian/geom/pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A/truth_H5977/trimesh/trimesh_0100_0050/trimesh';
% %     trimesh_file = fwd_shape;
% %     [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,nan,false,false);
% 
% %     elseif strcmp(fwd_shape,'prediction_as_truth')
% %         trimesh_file = './geom/TEST_population_sample_mean/trimesh_lung_0005';
% %         [tgeom,full_bdy] = read_2D_tri_mesh(trimesh_file,0.1,false,false);
% %         bspfile_mesh = "./geom/TEST_population_sample_mean/bspfile_torso.txt";
% %         bspfile_subject = "./geom/TEST_HLA-H11303_predicted/bspfile_torso.txt";
% %         [tgeom,~] = generate_subject_mesh(tgeom,full_bdy,bspfile_mesh,bspfile_subject); 
%     
% %     elseif strcmp(fwd_shape,'circle')
% %         mesh_shape = @(p)dcircle(p,0,0,150);                                           
% %         pfix = [];
% %         figure(101)
% %         [Nodes,~] = distmesh2d(mesh_shape,@huniform,h,[-1,-1;1,1],pfix); % plots
% %         full_bdy = Nodes(:,1).^2+Nodes(:,2).^2>0.999;                   % boundary nodes
% %         tgeom = delaunayTriangulation(Nodes);
% 
% %     else
% %         error
% %     end
% end
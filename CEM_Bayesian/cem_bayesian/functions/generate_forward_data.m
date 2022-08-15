function [data,sig_t_F,fwd_bdy,Le] = generate_forward_data(parameters,do_plots)

% addpath('distmesh','cem_bayesian')
rng(123)            

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param_m,param_e,param_d] = deal(parameters{:});
[mesh_file,seed_files] = deal(param_m{:});
[n_elec,n_patt,z_elec,w_elec,N,MeasPattern] = deal(param_e{:});
[nl,sig_0,conds] = deal(param_d{:});


%%%%%%%%%%%%%%%%%%%%%%%% Geometry & Meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [tgeom, full_bdy] = generate_mesh_fwd(fwd_shape);
[tgeom,full_bdy] = read_2D_tri_mesh(mesh_file,nan,false,false);
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

% Boundary
fwd_bdy = [x_bdy(ss_b),y_bdy(ss_b)];


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ground truth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Defined manually %%%%%%%%%%%%
%{obj_c,obj_w,obj_m};
% sig_t = mu_sig*ones(nn,1) + gaussian(obj_c{1},obj_w{1},obj_m{1},Nodes)+...
%                             gaussian(obj_c{2},obj_w{2},obj_m{2},Nodes);

%%%%%%%%%%% Defined from mesh seeds %%%%%%%%%%%
% Import lung shapes
xy_left = -readmatrix(seed_files{1});
xy_rght = -readmatrix(seed_files{2});
% Find nodes in lung shapes
poly_left = polyshape(xy_left(:,1),xy_left(:,2));
poly_rght = polybuffer(polyshape(xy_rght(:,1),xy_rght(:,2)), 1.);
xy_left_buff = polybuffer(poly_left, 0.1).Vertices;
xy_rght_buff = polybuffer(poly_rght, 0.1).Vertices;
lung_left = inpolygon(x,y,xy_left_buff(:,1),xy_left_buff(:,2));
lung_rght = inpolygon(x,y,xy_rght_buff(:,1),xy_rght_buff(:,2));
% Diseased states (conditions)
if ~isempty(conds)
    if any("LP"==conds), lung_left(y<mean(y)) = 0; end  % LEFT,  POSTERIOR
    if any("LA"==conds), lung_left(y>mean(y)) = 0; end  % LEFT,  ANTERIOR
    if any("RP"==conds), lung_rght(y<mean(y)) = 0; end  % RIGHT, POSTERIOR
    if any("RA"==conds), lung_rght(y>mean(y)) = 0; end  % RIGHT, ANTERIOR
end
% Define truth
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
params_e = {n_elec,in_elec,n_per_elec,z_elec,N};
params_p = {sig_true};
fem_params = {params_m,params_e,params_p};
[K,pK,SpK,i_K,M_B,C,D] = generate_fem_model(fem_params); 

% Solve
B = K+1/z_elec*M_B;
F = [zeros(nn,n_patt); N'*MeasPattern];
A = [B,C;C',D];
uu = A\F;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Meas_rest = sparse((1:n_elec-1),(nn+1:nn+n_elec-1),1);
data_t = MeasPattern*N*Meas_rest*uu;           % voltages at electrodes
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
if any(101==do_plots)
    plot_model('Forward model',101,0,tgeom,vars_bdy,vars_elec), hold on
    % plot(x(lung_left),y(lung_left),'or')
    % plot(x(lung_rght),y(lung_rght),'or')
    plot([xy_left(:,1);xy_left(1,1)],[xy_left(:,2);xy_left(1,2)],'k')
    plot([xy_rght(:,1);xy_rght(1,1)],[xy_rght(:,2);xy_rght(1,2)],'k')%, hold off
    fill([xy_left(:,1);xy_left(1,1)],[xy_left(:,2);xy_left(1,2)],[.7,.7,.7])
    fill([xy_rght(:,1);xy_rght(1,1)],[xy_rght(:,2);xy_rght(1,2)],[.7,.7,.7]), hold off
    kids = get(gca,'children');
    set(gca,'children',[kids(3:end);kids(1:2)])
end

% Plot the ground truth
if any(103==do_plots), plot_sigma('Ground truth',103,[1,2,1],tris,x,y,sig_true,0), end
if any(104==do_plots), plot_sigma('Ground truth',104,[3,2,1],tris,x,y,sig_true,0), end
if any(105==do_plots), plot_sigma('Ground truth',105,[3,3,1],tris,x,y,sig_true,0), end

end

function [] = run_eit_cem_bayes(eit_params,run_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                                                   %%%%%%%%%%%%%
%%%%%%%%%%%             COMPLETE ELECTRODE MODEL              %%%%%%%%%%%%%
%%%%%%%%%%%              BAYESIAN EIT INVERSION               %%%%%%%%%%%%%
%%%%%%%%%%%                                                   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root = '/hpc/mpag253/EIT/CEM_Bayesian/';
addpath('cem_bayesian/functions', ...
        'cem_bayesian/distmesh', ...
        'cem_bayesian/fem', ...
        'cem_bayesian/filexchange', ...
        'cem_bayesian/plot_functions')

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
rng(123) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNPACK
[params_m,params_e,params_d,params_p] = deal(eit_params{:});
[make_figs,save_figs,save_data,batch,num] = deal(run_params{:});
fprintf([repelem('=',1,60),'\n\nBATCH = ',batch,', NUMBER = ',num,'\n\n',repelem('=',1,60),'\n'])

% MESH
[pca_id,fwd_mesh,fwd_s1,fwd_s2,inv_mesh,inv_s1,inv_s2] = deal(params_m{:});
fwd_path = ['input/geom/',pca_id,'/',fwd_mesh];
inv_path = ['input/geom/',pca_id,'/',inv_mesh];
fwd_mesh_file = [fwd_path,'/trimesh/trimesh_',fwd_s1,'_',fwd_s2,'/trimesh.mat'];
inv_mesh_file = [inv_path,'/trimesh/trimesh_',inv_s1,'_',inv_s2,'/trimesh.mat'];
fwd_seed_files = {[fwd_path,'/mesh_seeds/mesh_seeds_lung_left_',fwd_s1,'.csv']; ...
                  [fwd_path,'/mesh_seeds/mesh_seeds_lung_right_',fwd_s1,'.csv']};


% ELECTRODES & PATTERNS
[n_elec,z_elec,w_elec] = deal(params_e{:});
n_patt = n_elec;                % # measurement patterns
% Basis functions for electrode potentials
N=-speye(n_elec);
N(:,1)=[];
N(1,:)=1;
% Current/measurement pattern (II)
MeasPattern=toeplitz([1;-1;zeros(n_elec-2,1)],[1,zeros(1,n_patt-2),-1]);

% DATA
[nl,conds] = deal(params_d{:});
sig_0 = 1;                      % (...)

% PRIORS
[prior_type,l_prcorr] = deal(params_p{:});
n_priors = 9;                   % # priors to plot

% RECONSTRUCTION
max_iter = 50;                  % max # of Gauss-Newton iterations


%%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nForward model:\n\n\t')

% Packing the input parameters
params_m = {fwd_mesh_file,fwd_seed_files};
params_e = {n_elec,n_patt,z_elec,w_elec,N,MeasPattern}; % <--- fix
params_d = {nl,sig_0,conds}; %,mu_sig}; % <--- fix
fwd_params = {params_m,params_e,params_d};

% Generating simulated data with forward solution
[data,sig_t_f,fwd_bdy,Le] = generate_forward_data(fwd_params,make_figs);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Inverse Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nInverse model:\n\n\t')

% Geometry & mesh
[tgeom,full_bdy] = read_2D_tri_mesh(inv_mesh_file,nan,true,false);
% Gamma_r = [];  % future work: include covariance of the node radii
nodes = tgeom.Points;
tris = tgeom.ConnectivityList;

% Generate the inverse model
params_m = {tgeom,nodes,tris,full_bdy};
params_e = {n_elec, w_elec};
inv_params = {params_m,params_e};
[model_vars] = generate_model(inv_params);
[vars_all,vars_bdy,vars_elec] = deal(model_vars{:});
[x,y,theta,nn] = deal(vars_all{:});
[x_bdy,y_bdy,theta_bdy,nn_bdy,bdy_indx,ss_b,bdy_elems] = deal(vars_bdy{:});
[elec_pts, in_elec, n_per_elec] = deal(vars_elec{:});

% Generate FEM model parameters
params_m = {nodes,tris,nn,bdy_indx,bdy_elems};
params_e = {n_elec,in_elec,n_per_elec,z_elec,N};
params_d = {sig_0};
fem_params = {params_m,params_e,params_d};
[~,pK,SpK,i_K,M_B,C,D] = generate_fem_model(fem_params);

% Forcing function
F = [zeros(nn,n_patt);N'*MeasPattern];

% (...)
B_imp = 1/z_elec*M_B;

% Plot inverse model
if any(102==make_figs)
    plot_model('Inverse model',102,0,tgeom,vars_bdy,vars_elec)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create priors, mean shape info, and regularisation matrix (L_sig)

% Priors derived from statistical shape model
params_m = {nn,tris,x,y,theta};
params_p = {prior_type,n_priors,l_prcorr};
params_f = {root,pca_id,inv_mesh_file};
prior_params = {params_m,params_p,params_f};
[mu_sig,L_sig,G_sig] = generate_priors(prior_params,false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate solution
sig_new = mu_sig; %*ones(nn,1);
iter = 1;
norm_grad_init = 1;
Meas_rest = sparse((1:n_elec-1),(nn+1:nn+n_elec-1),1);
data_v = data(:);
tol = 1e-7; % 1e-7; % 
converged = false;

% Plot the initial reconstruction
clims1 = [0,1]; %[min(sig_t) max(sig_t)];
if any(103==make_figs)
    plot_sigma('Initial estimate',103,[1,2,2],tris,x,y,sig_new,clims1)
end

% Iterate solution
fprintf('\nReconstruction:\n\n')
while (iter<=max_iter) && ~converged
    
    K = sparse(i_K(:,1),i_K(:,2),pK*exp(sig_new),nn,nn);
    B = K+B_imp;
    A = [B,C;C',D];
    uu_new = A\F;
    d_L = MeasPattern*N*(Meas_rest/A); 

    S = sparse(n_patt*n_elec,nn);
    for jj = 1:n_patt
        S(1+(jj-1)*n_elec:jj*n_elec,:) = ...
            -d_L*[reshape(SpK*uu_new(1:nn,jj),nn,nn);zeros(n_elec-1,nn)];
    end

    J_total = [Le*S*diag(exp(sig_new));L_sig];
    b_total = [Le*(reshape(MeasPattern*N*Meas_rest*uu_new,n_patt*n_elec,1)-data_v);...
               L_sig*(sig_new-mu_sig)];
    dir_k = -J_total\b_total;

    % Cost Function
    P = @(kappa) norm(Le*(reshape(MeasPattern*N*Meas_rest*([sparse(i_K(:,1),i_K(:,2),pK*...
        exp(sig_new+kappa*dir_k))+B_imp,C;C',D]\F),n_patt*n_elec,1)-data_v))^2+...
        norm(L_sig*(sig_new+kappa*dir_k-mu_sig))^2;

    % Line Search
    % len=1;    % to not do line search
    len = fminsearch(P,1);
    sig_new = sig_new+len*dir_k;
    norm_grad = norm(J_total'*b_total);

    % Update iterations
    fprintf('\tIteration %2d: norm grad = %s\n',iter,num2str(norm_grad))
    if iter == 1, norm_grad_init = norm_grad; end
    iter = iter+1;
    if norm_grad/norm_grad_init<tol, converged = true; fprintf('\tConverged.\n'), end
    if iter > max_iter, fprintf('\tFailed to converge.\n'); end
    
    % Plot latest reconstruction
    if any(103==make_figs)
        if converged, pt = 'Reconstruction'; else, pt = 'Latest iteration'; end
        plot_sigma(pt,103,[1,2,2],tris,x,y,sig_new,clims1)
        pause(.5)
    end
end

% Position figure
scrnwh = get(0,'screensize'); fig_w = 1200; fig_h = 600;
set(gcf,'Position',[0 scrnwh(4)-fig_h fig_w fig_h])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POSTERIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nPosterior:\n')

%  Posterior calculations
G_post = inv(J_total'*J_total);               % posterior covariance
S_post = diag(G_post).^(.5);                  % standard deviations
L_post = chol(G_post,'lower');                % (...)

% Samples from the posterior distrubution
n_post_samps = 1e4;
sig_samps = sig_new+L_post*randn(nn,n_post_samps);

% Ground truth evaluated on inverse solution mesh
% sig_t = sig_t_f(x,y);  % absolute
sig_t = get_relative_ground_truth(sig_t_f,x,y,x_bdy,y_bdy,fwd_bdy);  %relative

% Relative error of solution
% rel_err = norm(sig_new-sig_t)/norm(sig_t);
rel_err = norm(get_truncated_sigma(sig_new)-sig_t)/norm(sig_t);

%  Mahalanobis distance of solution
m_dist = get_mahalanobis_distance(x,y,tris,sig_new,sig_t,G_post);

%  Total variation
[tv_mean,tv_stdv] = get_total_variation(nodes,tris,x,y,sig_new,sig_samps);
% alternative
tv_alt = trace(G_post);
% overall_correlation = det(G_post);

% Display
fprintf("\n\tPrior corr. len.:   %.2e",l_prcorr)
fprintf("\n\tRelative error:     %.2e",rel_err)
fprintf("\n\tMahalanobis dist.:  %.2e",m_dist)
fprintf("\n\tTotal variation:    %.2e \x00B1 %.2e\n",tv_mean,tv_stdv)
fprintf("\n\t(Total variation ALT):    %.2e\n",tv_alt)
fprintf("\n\t(overall correlation):    %.2e\n",det(G_post))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Summary plots
if any(104==make_figs)
    plot_sigma('Prior Mean',104,[3,2,3],tris,x,y,mu_sig,clims1)
    plot_sigma('Prior SD',104,[3,2,5],tris,x,y,diag(G_sig).^.5,clims1)
    plot_sigma('MAP Estimate',104,[3,2,4],tris,x,y,sig_new,clims1)
    plot_sigma('Posterior SD',104,[3,2,6],tris,x,y,S_post,clims1)
    plot_metrics(104,[3,2,2],l_prcorr,m_dist,tv_mean,tv_stdv)
end

% Solution slice
if any(105==make_figs)
    plot_sigma('MAP Estimate',105,[3,3,2],tris,x,y,sig_new,clims1)
    scrnwh = get(0,'screensize'); fig_w = 1000; fig_h = 750;
    set(gcf,'Position',[0 scrnwh(4)-fig_h fig_w fig_h])
    plot_metrics(105,[3,3,3],rel_err,m_dist,tv_mean,tv_stdv)
    plot_add_section('A',105,[3,3,2],x,-30)
    plot_add_section('B',105,[3,3,2],x, 30)
    ax_A = plot_section('Section A',105,[9,2,7:2:15], x,y,-30,sig_t,sig_new,S_post,sig_samps(:,1:3));
    ax_B = plot_section('Section B',105,[9,2,8:2:16], x,y, 30,sig_t,sig_new,S_post,sig_samps(:,1:3));
    linkaxes([ax_A ax_B], 'y')
    % font size
    kids=get(gcf,'children');
    for i = 1:length(kids), set(kids(i),'fontsize', 14); end
    % legend
    kids=get(ax_A, 'Children');
    lgd = legend([kids(2) kids(1) kids(6) kids(5)], ...
                 {'True {} ','MAP {} ','${99\%}$ CI {} ','Samples'}, ...
                 'Orientation','horizontal');
    lgd_pos = lgd.Position;
    lgd_pos(1) = (1-lgd_pos(3))/2;
    lgd_pos(2) = 0.075;
%     lgd.Position(lgd_pos);
    set(lgd,'Position',lgd_pos)

end

% Posterior samples
if any(106==make_figs)
'eit_animated.gif'    fname = 'eit_animated.gif'
    plot_post_samps('Posterior Samples',106,0,tris,x,y,sig_samps(:,1:50),clims1,fname)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nOutput:\n\n')
if isempty(save_figs) && ~save_data, fprintf('\tNone.\n'); end

% Data
if save_data
    fname = ['output/metrics/metrics_',batch,'.xlsx'];
    fhead = {'Batch','Num','Converged','PCL','MD','TV','TV (sd)'};
    fdata = {batch,num,converged,l_prcorr,m_dist,tv_mean,tv_stdv};
    if ~isfile(fname),  writecell(fhead,fname); end
    writecell(fdata,fname,'WriteMode','append');
    fprintf(['\tSaved "',fname,'"\n'])
end

% Save plots
for fg = save_figs
    save_dir = ['output/figures/',batch,'/',batch,'-',num2str(num,'%03d'),'/'];
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    fname = [save_dir,'figure_',num2str(fg),'.png'];
    figure(fg), exportgraphics(gcf,fname,'Resolution',150)
    fprintf(['\tSaved "',fname,'"\n'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
end

% function [m_dist] = get_mahalanobis_distance(x,y,tris,sig_new,sig_t,G_post)
%     % weight mahalanobis distance using element areas
%     md_wts = zeros(length(x),1);
%     elem_areas = zeros(size(tris,1),1);
%     for e = 1:size(tris,1)
%         elem_areas(e) = polyarea(x(tris(e,:)),y(tris(e,:)));
%     end
%     elem_areas = elem_areas/sum(elem_areas);
%     for n = 1:length(x)
%         nod_elems = [find(tris(:,1)==n);find(tris(:,2)==n);find(tris(:,3)==n)];
%         md_wts(n) = sum(elem_areas(nod_elems))/3;
%     end
%     % figure(999), trisurf(tris,x,y,md_wts)
%     % size(md_wts)
%     % size((sig_t-sig_new)')
%     wtd_diff = md_wts.*(sig_t-sig_new);
%     % size(wtd_diff)
%     m_dist = (wtd_diff'/G_post*wtd_diff)^0.5;
%     % m_dist = ((sig_t-sig_new)/(G_post)*(sig_t-sig_new))).^0.5;
%     % p-value: d^2 follows the chi-squared distribution with n degrees of freedom, 
%     % where n is the number of dimensions of the normal distribution. n = ?
% end

% function [sig_trunc] = get_truncated_sigma(sig)
%     sig_trunc = sig;
%     sig_trunc(sig_trunc>1) = 1;
%     sig_trunc(sig_trunc<0) = 0;
% end



function [] = run_eit_cem_bayes(eit_params)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%                                                   %%%%%%%%%%%%%
    %%%%%%%%%%%             COMPLETE ELECTRODE MODEL              %%%%%%%%%%%%%
    %%%%%%%%%%%              BAYESIAN EIT INVERSION               %%%%%%%%%%%%%
    %%%%%%%%%%%                                                   %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(groot,'defaulttextinterpreter','latex');  
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    set(groot,'defaultLegendInterpreter','latex');

    root = '/hpc/mpag253/EIT/CEM_Bayesian/';
    addpath('cem_bayesian/functions', ...
            'cem_bayesian/distmesh', ...
            'cem_bayesian/fem', ...
            'cem_bayesian/filexchange', ...
            'cem_bayesian/plot_functions')
    rng(123) 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % UNPACK
    [params_m,params_e,params_d,params_p] = deal(eit_params{:});
    
    % MESH
    [pca_id,fwd_mesh,inv_mesh] = deal(params_m{:});
    fwd_mesh_file = strcat(root,'geom/',pca_id,'/',fwd_mesh,'/trimesh');
    inv_mesh_file = strcat(root,'geom/',pca_id,'/',inv_mesh,'/trimesh');
    
    % ELECTRODES & PATTERNS
    [n_elec,n_patt,z_elec,w_elec] = deal(params_e{:});
    % (...)
    N_=-speye(n_elec);
    N_(:,1)=[];
    N_(1,:)=1;
    % Current/measurement pattern (II)
    MeasPattern=toeplitz([1;-1;zeros(n_elec-2,1)],[1,zeros(1,n_patt-2),-1]);

    % DATA
    [nl,sig_0] = deal(params_d{:});

    % PRIORS
    [~,~,l_prcorr] = deal(params_p{:});

    % RECONSTRUCTION
    max_iter = 50;                  % max # of Gauss-Newton iterations
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nForward model:\n\n\t')
    
    % Packing the input parameters
    params_m = {fwd_mesh_file};
    params_e = {n_elec,n_patt,z_elec,w_elec,N_,MeasPattern}; % <--- fix
    params_d = {nl,sig_0}; %,mu_sig}; % <--- fix
    fwd_params = {params_m,params_e,params_d};
    
    % Generating simulated data with forward solution
    [data,sig_t_f,Le] = generate_forward_data(fwd_params);
    
    
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
    params_e = {n_elec,in_elec,n_per_elec,z_elec,N_};
    params_d = {sig_0};
    fem_params = {params_m,params_e,params_d};
    [~,pK,SpK,i_K,M_B,C,D] = generate_fem_model(fem_params);

    % Forcing function
    F = [zeros(nn,n_patt);N_'*MeasPattern];
    
    % (...)
    B_imp = 1/z_elec*M_B;
    
    % Plot inverse model
    plot_model('Inverse model',124,0,tgeom,vars_bdy,vars_elec)
    % % set(gcf,'Position',[f4w scrnwh(4)-f3h f3w f3h])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create priors, mean shape info, and regularisation matrix (L_sig)
    
    % Priors derived from statistical shape model
    params_m = {nn,tris,x,y,theta};
    %params_p = {prior_type,n_priors,l_prcorr};
    params_f = {root,pca_id,inv_mesh};
    prior_params = {params_m,params_p,params_f};
    [mu_sig, L_sig, G_sig] = generate_priors(prior_params,false);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initiate solution
    sig_new = mu_sig; %*ones(nn,1);
    iter = 1;  %0;
    norm_grad = 1;
    norm_grad_init = 1;
    Meas_rest = sparse((1:n_elec-1),(nn+1:nn+n_elec-1),1);
    data_v = data(:);
    tol = 1e-7; % 1e-7; % 
    
    % Plot the initial reconstruction
    clims1 = [0,1]; %[min(sig_t) max(sig_t)];
    plot_sigma('Initial estimate',103,[1,2,2],tris,x,y,sig_new,clims1)
    
    % Iterate solution
    fprintf('\nReconstruction:\n\n')
    while iter <= max_iter && norm_grad/norm_grad_init>tol
        
        K = sparse(i_K(:,1),i_K(:,2),pK*exp(sig_new),nn,nn);
        B = K+B_imp;
        A = [B,C;C',D];
        uu_new = A\F;
        d_L = MeasPattern*N_*(Meas_rest/A); 
    
        S = sparse(n_patt*n_elec,nn);
        for jj = 1:n_patt
            S(1+(jj-1)*n_elec:jj*n_elec,:) = ...
                -d_L*[reshape(SpK*uu_new(1:nn,jj),nn,nn);zeros(n_elec-1,nn)];
        end
    
        J_total = [Le*S*diag(exp(sig_new));L_sig];
        b_total = [Le*(reshape(MeasPattern*N_*Meas_rest*uu_new,n_patt*n_elec,1)-data_v);...
                   L_sig*(sig_new-mu_sig)];
        dir_k = -J_total\b_total;
    
        % Cost Function
        P = @(kappa) norm(Le*(reshape(MeasPattern*N_*Meas_rest*([sparse(i_K(:,1),i_K(:,2),pK*...
            exp(sig_new+kappa*dir_k))+B_imp,C;C',D]\F),n_patt*n_elec,1)-data_v))^2+...
            norm(L_sig*(sig_new+kappa*dir_k-mu_sig))^2;
    
        % Line Search
        % len=1;    % to not do line search
        len = fminsearch(P,1);
        sig_new = sig_new+len*dir_k;
        norm_grad = norm(J_total'*b_total);
    
        % Update iterations
        fprintf('\titeration %2d: norm grad = %s\n',iter,num2str(norm_grad))
        if iter == 1, norm_grad_init = norm_grad; end
        iter = iter+1;
        
        % Plot latest reconstruction
        plot_sigma('Latest iteration',103,[1,2,2],tris,x,y,sig_new,clims1)
        pause(.5)
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POSTERIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPosterior:\n')
    
    %  Posterior calculations
    G_post = inv(J_total'*J_total);               % posterior covariance
    S_post = diag(G_post).^(.5);                % standard deviations
    L_post = chol(G_post,'lower');                % (...)
    
    % Samples from the posterior distrubution
    n_post_samps = 1e2;
    sig_samps = sig_new+L_post*randn(nn,n_post_samps);
    
    % Ground truth evaluated on inverse solution mesh
    sig_t = sig_t_f(x,y);
    
    %  Mahalanobis distance of solution
    m_dist = (((sig_t-sig_new)')/(G_post)*(sig_t-sig_new)).^0.5;
    % p-value: d^2 follows the chi-squared distribution with n degrees of freedom, 
    % where n is the number of dimensions of the normal distribution.
    %  n = ?
    
    %  Total variation
    [tv_mean,tv_stdv] = get_total_variation(nodes,tris,x,y,sig_new,sig_samps);
    
    % Display
    metrics = sprintf("\n\tPrior corr. len.:   %.2e\n" + ...
                       "\n\tMahalanobis dist.:  %.2e\n" + ...
                       "\n\tTotal variation:    %.2e \x00B1 %.2e\n\n", ...
                       l_prcorr,m_dist,tv_mean,tv_stdv);
    fprintf(metrics)
    % fprintf("\n\tPrior corr. len.:   %.2e\n",l_prcorr)
    % fprintf("\n\tMahalanobis dist.:  %.4f\n",m_dist)
    % fprintf("\n\tTotal variation:    %.4f \x00B1 %.4f\n\n",tv_mean,tv_stdv)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Setup
    % scrnwh = get(0,'screensize');
    % f5w = 900; f5h = 400;
    clims2 = [min(min(sig_samps)) max(max(sig_samps))];
    
    % Summary plots
    % plot_sigma('Ground Truth',104,[3,2,[1,2]],)
    plot_sigma('Prior Mean',104,[3,2,3],tris,x,y,mu_sig,clims1)
    plot_sigma('Prior SD',104,[3,2,5],tris,x,y,diag(G_sig).^.5,clims1)
    plot_sigma('MAP Estimate',104,[3,2,4],tris,x,y,sig_new,clims1)
    plot_sigma('Posterior SD',104,[3,2,6],tris,x,y,S_post,clims1)
    metrics = {'Metrics','', ...
               ['PCL:{} ',num2str(l_prcorr,'%.2e')],'', ...
               ['MD: {} ',num2str(m_dist,  '%.2e')],'', ...
               ['TV: {} ',num2str(tv_mean, '%.2e'),' PM ',num2str(tv_stdv,'%.2e')]};
    metrics{1} = ['\underline{\bf ', metrics{1}, ' \rm}'];
    metrics{7} = strrep(metrics{7},'PM','$\pm$');
    figure(104), subplot(3,2,2)
    text(.5,.5,metrics, ...
        'fontsize',12, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle')
    set(gca,'visible','off')
    drawnow
    
    % Solution slice
    x_s = linspace(min(x),max(x),1e3);    % section x-coordinates 
    y_s1 = -30*ones(size(x_s));           % section 1 y-coordinates
    y_s2 =  30*ones(size(x_s));           % section 2 y-coordinates
    plot_sigma('MAP Estimate',105,[3,2,[1,2]],tris,x,y,sig_new,clims1)
    plot_add_section('A',105,[3,2,[1,2]],x_s,y_s1)
    plot_add_section('B',105,[3,2,[1,2]],x_s,y_s2)
    s1 = get_sections(x,y,x_s,y_s1,sig_t,sig_new,S_post,sig_samps(:,1:3)); % <--- combine with plot_section
    s2 = get_sections(x,y,x_s,y_s2,sig_t,sig_new,S_post,sig_samps(:,1:3));
    plot_section('Section A',105,[3,2,[3,5]],x_s,s1)
    plot_section('Section B',105,[3,2,[4,6]],x_s,s2)
    % set(gcf,'Position',[scrnwh(3)-f5w scrnwh(4)-f5h f5w f5h])
    
    % % Posterior samples
    % plot_post_samps('Posterior Samples',1,0,tris,x,y,sig_samps,clims1,'eit_animated.gif')
    
    % % Level set
    % sig_new_ls = get_level_set(sig_new);
    % plot_sigma('Prior Mean',42,0,tris,x,y,sig_new_ls,0)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     function [data_LS] = get_level_set(data)
%         data_LS=0*data;
%         data_LS(data>0.166) = 0.333;
%         data_LS(data>0.500) = 0.666;
%         data_LS(data>0.833) = 1.00;
%     end
    
%     function [tv_mean,tv_stdv] = get_total_variation(nodes,tris,x,y,sig_new,sigma_samples)
%         %  first calculate mesh element areas
%         tri_areas = zeros([length(tris) 1]);
%         for i = 1:length(tri_areas)
%             x1 = nodes(tris(i,1),1); y1 = nodes(tris(i,1),2);
%             x2 = nodes(tris(i,2),1); y2 = nodes(tris(i,2),2);
%             x3 = nodes(tris(i,3),1); y3 = nodes(tris(i,3),2);
%             tri_areas(i) = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
%         end
%         %  then gradient and total variance at each element (relative to mean) 
%         %  for each posterior sample
%         tot_vars = zeros([size(sigma_samples,2) 1]);
%         for ii=1:size(sigma_samples,2)
%             [dFx,dFy] = trigradient(tris,x,y,sigma_samples(:,ii)-sig_new,'face'); % at faces!
%             tot_vars(ii) = sum((abs(dFx)+abs(dFy)).*tri_areas)/sum(tri_areas);
%         end
%         tv_mean = mean(tot_vars);
%         tv_stdv = std(tot_vars);
%     end




end
function [mu_sig,L_sig,G_sig] = generate_priors(prior_params,do_plot)
%generate_priors Creates priors and returns regularisation matrix (L)
%   ..

% Unpack
[params_m,params_p,params_f] = deal(prior_params{:});
[nn,Tris,x,y,theta] = deal(params_m{:});
[prior_id,npriors,prior_corr_length] = deal(params_p{:});
[path,pca_id,inv_mesh] = deal(params_f{:});

% Load probability information
mu_sig_file = strcat(path,'geom/',pca_id,'/',inv_mesh,'/nodal_pf_mean_s-10000.csv');
c_lung_file = strcat(path,'geom/',pca_id,'/',inv_mesh,'/nodal_pf_condcov_s-10000.csv');
mu_sig_orig = readmatrix(mu_sig_file);
c_lung_orig = readmatrix(c_lung_file);

% % Plot sample slice of the probability imformation
% x_samples = linspace(min(x),max(x),1e3); y_samples = 0*x_samples;
% F = scatteredInterpolant(x,y,mu_sig_orig); z1 = F(x_samples,y_samples);
% F = scatteredInterpolant(x,y,mu_sig_orig+diag(c_lung_orig)); z2 = F(x_samples,y_samples);
% F = scatteredInterpolant(x,y,mu_sig_orig-diag(c_lung_orig)); z3 = F(x_samples,y_samples);
% figure(999), hold on
% plot(x_samples,z1,'k-')
% plot(x_samples,z2,'k--')
% plot(x_samples,z3,'k--')

% For separate lung probabilities i.e. no covariance between lungs
mu_sig = mu_sig_orig;
c_lung = c_lung_orig;

% % For combined lung probabilities: modify the variance/covariance to
% % include 0
% mu_sig = mu_sig_orig./2;
% c_lung = 0.5*(c_lung_orig + mu_sig*mu_sig');
% % issymmetric(c_lung)
% % min(eig(c_lung))
% % Plot sample slice of modified probability functions
% x_samples = linspace(min(x),max(x),1e3); y_samples = 0*x_samples;
% F = scatteredInterpolant(x,y,mu_sig); z1 = F(x_samples,y_samples);
% F = scatteredInterpolant(x,y,mu_sig+diag(c_lung).^.5); z2 = F(x_samples,y_samples);
% F = scatteredInterpolant(x,y,mu_sig-diag(c_lung).^.5); z3 = F(x_samples,y_samples);
% figure(998), hold on
% plot(x_samples,z1,'k-')
% plot(x_samples,z2,'k--')
% plot(x_samples,z3,'k--')
% hold off
% stop

% Gaussians for each node
dnsq = (x-x').^2+(y-y').^2;
G_sig_1 = exp(-1/prior_corr_length^2*dnsq);
% dnsq = abs(x-x')+abs(y-y');
% G_sig_1 = exp(-1/prior_corr_length*dnsq);

% Prior factor
prior_factor = 1.;


% Create the priors

% p0: No shape information, pure gaussian correlations
if strcmp(prior_id,'p0')
    G_sig = prior_factor*G_sig_1;
    mu_sig=0.5*(0*mu_sig+1);

% p1: Mean shape information, pure gaussian correlations
elseif strcmp(prior_id,'p1')
    G_sig = prior_factor*G_sig_1;

% p2: Gaussian correlations with shape information
elseif strcmp(prior_id,'p2')
    G_sig_2 = c_lung;
    G_sig = prior_factor*G_sig_1.*G_sig_2;

% p3: Gaussian correlations with shape information, penalising non-lung    
elseif strcmp(prior_id,'p3')
    G_sig_2 = c_lung;
    G_sig_3 = (mu_sig_orig*mu_sig_orig').^0.5;
    G_sig = prior_factor*G_sig_1.*G_sig_2.*G_sig_3;

else
    error(strcat(['Invalid prior identifier: ' prior_id]))
end

% Add small values to diag
G_sig = G_sig + 1e-6*eye(nn);

% Regularisation matrix
L_sig = chol(inv(G_sig));

% Test compatibility of input data
if size(mu_sig,1) ~= size(L_sig,1)
    error('Incompatible data for priors. Covariance based on the wrong mesh?')
end

% % Plot the priors
scrnwh = get(0,'screensize');
if do_plot
    rd = randn(nn,npriors);             % random for each node and prior
    draw = mu_sig+L_sig\rd;             % to visualise the priors
    subpltdim = ceil(sqrt(npriors));    % dimension of subplot to show priors
    for ii = 1:npriors
        figure(102); subplot(subpltdim,subpltdim,ii)
        f2w = 600; f2h = 500; set(gcf,'Position',[scrnwh(3)-f2w scrnwh(4)-f2h f2w f2h])
        trisurf(Tris,x,y,draw(:,ii))
        view(2), shading interp, daspect([1,1,1]), axis off
        cb = colorbar; cb.TickLabelInterpreter = 'latex';
        caxis([min(draw(:)),max(draw(:))]), set(gca,'fontsize',10)
        title(strcat('prior sample',{' '},int2str(ii)),'interpreter','latex')
        drawnow
    end
end

% % Plot G_sig_2
% for i = [100 200 300 400 500]%[1485,194,157,1077] % y(find((x>50) & (x<51)))
%     figure(2e4+i), trisurf(Tris,x,y,G_sig_2(:,i))
%     view(2), shading interp, axis equal, axis off
% %     caxis([0,1])
% end
% % % Plot G_sig
% % for i = [100 200 300 400 500]%[1485,194,157,1077] % y(find((x>50) & (x<51)))
% %     figure(1e4+i), trisurf(Tris,x,y,G_sig(:,i))
% %     view(2), shading interp, axis equal, axis off
% % %     caxis([0,1])
% % end
% stop

end
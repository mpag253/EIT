function [mu_sig,L_sig] = generate_priors(nn,Tris,x,y,theta,npriors,prior_corr_length,do_plot)
%generate_priors Creates priors and returns regularisation matrix (L)
%   ..

% % % Probability information
mu_sig_orig = readmatrix('./geom/TEST_HLA-H11303_predicted/TEST_lung_nodes_mean.csv');
% mu_sig = 0*mu_sig
size(mu_sig_orig)
min(mu_sig_orig),max(mu_sig_orig) % PRINTS!
c_lung_orig = readmatrix('./geom/TEST_HLA-H11303_predicted/TEST_lung_nodes_cov.csv');
% % Plot
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
% % Plot
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

% Create the priors

% Gaussians for each node
dnsq = (x-x').^2+(y-y').^2;
G_sig_1 = exp(-1/prior_corr_length^2*dnsq);
% dnsq = abs(x-x')+abs(y-y');
% G_sig_1 = exp(-1/prior_corr_length*dnsq);

% Shape information
G_sig_2 = c_lung;
% min(min(G_sig_2)),max(max(G_sig_2)),min(diag(G_sig_2)),max(diag(G_sig_2)) % PRINTS!

% Lung factor
G_sig_3 = (mu_sig_orig*mu_sig_orig').^0.5;
% G_sig_3 = (mu_sig_orig + mu_sig_orig')./2; %DOESN'T WORK!
% min(min(G_sig_3)),max(max(G_sig_3)),min(diag(G_sig_3)),max(diag(G_sig_3)) % PRINTS!

% Prior factor
prior_factor = 1.;

% Combined
G_sig = prior_factor*G_sig_1;                  % no shape info (original)
% G_sig = prior_factor*G_sig_1.*G_sig_2;         % proximity + shape
% G_sig = prior_factor*G_sig_1.*G_sig_2.*G_sig_3;   % proximity + shape + lung factor

% Add small values to diag
G_sig = G_sig + 1e-6*eye(nn);


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


% Regularisation matrix (?)
% Upper triangular of (?) sigma (Cholesky factorisation)
L_sig = chol(inv(G_sig));


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

end
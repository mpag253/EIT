function [ax] = plot_section(pt,fg,sp,x,y,y_s,sig_t,sig_new,S_post,sig_samps)

    % Define section
    x_s = linspace(min(x),max(x),1e3);    % section x-coordinates 
    y_s = y_s*ones(size(x_s));            % section 1 y-coordinates

    % Truth
    F = scatteredInterpolant(x,y,sig_t);
    z_t = F(x_s,y_s);
    
    % MAP estimate
    F = scatteredInterpolant(x,y,sig_new);         
    z_s = F(x_s,y_s);
    
    % Posterior standard devaition
    F = scatteredInterpolant(x,y,S_post);          
    z_s_sd = F(x_s,y_s);
    
    % Posterior samples
    n_samps = size(sig_samps,2);
    z_s_samps = cell(n_samps);
    for i = 1:n_samps
        F = scatteredInterpolant(x,y,sig_samps(:,i));
        z_s_samps{i} = F(x_s,y_s);
    end
    
    % Confidence interval
    z_s_ci = {z_s+2.576*z_s_sd, z_s-2.576*z_s_sd};
    

    %%

    % Figure setup
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    hold on

    % Confidence interval
    fill([x_s,       fliplr(x_s)], ...
         [z_s_ci{1}, fliplr(z_s_ci{2})], ...
         .8*[1 1 1], 'LineStyle','none', 'FaceAlpha',0.5);
    % plot(x_s,z_s_ci(1),'--','LineWidth',2,'Color',1*[1 1 1]);
    % plot(x_s,z_s_ci(2),'--','LineWidth',2,'Color',1*[1 1 1]);
    
    % Posterior samples
    for i = 1:length(z_s_samps)
        plot(x_s,z_s_samps{i},':','LineWidth',2,'Color',.5*[1 1 1]);
    end

    % Truth & MAP estimate
    plot(x_s,z_t,'r-','LineWidth',2);
    plot(x_s,z_s,'-','LineWidth',2,'Color',.0*[1 1 1]);

    % Format plot
    xlim([min(x_s) max(x_s)])
    xlabel('x [mm]')
    ylabel('${\sigma}$')
    ytickformat('%.1f')
    pt = ['\bf ', pt, ' \rm'];
    title(pt,'interpreter','latex')
    plots=get(gca, 'Children');
%     legend([plots(2) plots(1) plots(6) plots(5)], ...
%            {'True','MAP','${99\%}$ CI','Samples'}, ...
%            'Location','SouthOutside', ...
%            'Orientation','horizontal')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    set(gca, 'layer', 'top')
    set(gca, 'fontsize', 14)
    hold off

    ax = gca;
end
function [] = plot_section(pt,fg,sp,x_s,section_data)

    % Unpack section data
    [z_t,z_s,z_s_samps,z_s_ci] = deal(section_data{:});
    
    % Setup
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
    xlabel('x')
    ylabel('${\sigma}$')
    title(pt,'interpreter','latex')
    plots=get(gca, 'Children');
    legend([plots(2) plots(1) plots(6) plots(5)], ...
           {'True','MAP','${99\%}$ CI','Samples'}, ...
           'Location','SouthOutside', ...
           'Orientation','horizontal')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    set(gca, 'layer', 'top')
    hold off
end
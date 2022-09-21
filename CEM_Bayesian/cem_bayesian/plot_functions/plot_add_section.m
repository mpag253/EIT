function [] = plot_add_section(label,fg,sp,x,y_s)
    x_s = linspace(min(x),max(x),1e3);    % section x-coordinates 
    y_s = y_s*ones(size(x_s));            % section 1 y-coordinates
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    hold on
    plot3(x_s,y_s,1e6*ones(length(x_s)),'k-')
    text(max(x_s),max(y_s),1e6,label, ...
         'FontSize',14, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle')
    hold off
end
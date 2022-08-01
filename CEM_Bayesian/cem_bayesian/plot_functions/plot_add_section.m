function [] = plot_add_section(label,fg,sp,x_s,y_s)
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    hold on
    plot3(x_s,y_s,1e6*ones(length(x_s)),'k-')
    text(max(x_s),max(y_s),1e6,label, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle')
    hold off
end
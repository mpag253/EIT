function [] = plot_sigma(plot_title,fg,sp,tris,x,y,sigma,clims)
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    trisurf(tris,x,y,sigma)
    axis off
    daspect([1,1,1])
    if size(clims,2)>1, caxis(clims), end
    % xlim([-250 250]), ylim([-200 200])
    view(2)
    shading interp
    set(gca,'fontsize',10)
    title(plot_title,'interpreter','latex')
    drawnow
end
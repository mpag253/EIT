function [] = plot_tau(pt,fg,sp,tris,x,y,tau,clims)
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    colormap(coolwarm)
    trisurf(tris,x,y,tau)
    axis off
    daspect([1,1,1])
    if size(clims,2)>1, caxis(clims), end
    % xlim([-250 250]), ylim([-200 200])
    view(2)
    shading interp
    set(gca,'fontsize',14)
    pt = ['\bf ', pt, ' \rm'];
    title(pt) %,'interpreter','latex')
    drawnow
end
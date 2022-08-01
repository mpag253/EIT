function [] = plot_model(pt,fg,sp,tgeom,vars_bdy,vars_elec)

    % Unpack
    [x_bdy,y_bdy,~,~,~,ss_b,~] = deal(vars_bdy{:});
    [elec_pts, in_elec, ~] = deal(vars_elec{:});

    % Setup
    figure(fg)
    if sp, subplot(sp(1),sp(2),sp(3:end)), end
    hold on

    % Plot mesh
    triplot(tgeom, 'Color', [.5 .5 .5])

    % Plot boundary
    plot(x_bdy(ss_b),y_bdy(ss_b),'-k.')

    % Plot electrode nodes
    in_any_elec = logical(sum(in_elec,2));
    plot(x_bdy(in_any_elec),y_bdy(in_any_elec),'og','MarkerSize',4,'LineWidth',2)

    % Plot electrode points
    plot(elec_pts(:,1),elec_pts(:,2),'xb','MarkerSize',12,'LineWidth',3)

    % Plot electrode labels
    spacing_factor = 1.1;
    for i = 1:length(elec_pts)
        text(1.2*elec_pts(i,1),spacing_factor*elec_pts(i,2),num2str(i),...
             'Color','b','FontSize',10)
    end

    % Format
    axis off, daspect([1,1,1])
    view(2), %set(gca,'fontsize',10)
    title(pt,'interpreter','latex')
    hold off
end
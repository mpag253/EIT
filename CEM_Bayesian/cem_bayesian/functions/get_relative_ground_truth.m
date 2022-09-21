function [sig_t] = get_relative_ground_truth(sig_t_f,x,y,x_bdy,y_bdy,fwd_bdy)
    % get polar of forward shape
    [t_fwd_bdy,r_fwd_bdy] = cart2pol(fwd_bdy(:,1),fwd_bdy(:,2));
    [t_fwd_bdy,ss_fwd] = sort(t_fwd_bdy);
    r_fwd_bdy = r_fwd_bdy(ss_fwd);
    % get polar of inverse boundary
    [t_inv_bdy,r_inv_bdy] = cart2pol(x_bdy,y_bdy);
    [t_inv_bdy,ss_inv] = sort(t_inv_bdy);
    r_inv_bdy = r_inv_bdy(ss_inv);
    % get polar of inverse nodes
    [t_nodes,r_nodes] = cart2pol(x,y);
    % boundary scale functions
    r_nodes_fwd_bdy = interp1(t_fwd_bdy,r_fwd_bdy,t_nodes,'spline');
    r_nodes_inv_bdy = interp1(t_inv_bdy,r_inv_bdy,t_nodes,'spline');
    % scale mesh
    r_nodes_scaled = (r_nodes./r_nodes_inv_bdy).*r_nodes_fwd_bdy;
    % revert to cartesian
    [x_fwd_map,y_fwd_map] = pol2cart(t_nodes,r_nodes_scaled);
    % evaluate truth
    sig_t = sig_t_f(x_fwd_map,y_fwd_map);
    % % plot
    % figure(999), hold on
    % % polarplot(t_fwd_bdy,r_fwd_bdy,'ro')
    % % polarplot(t_nodes,r_nodes,'b.')
    % % % polarplot(t_nodes,r_nodes_fwd_bdy,'rx')
    % % % polarplot(t_nodes,r_nodes_inv_bdy,'bx')
    % % polarplot(t_nodes,r_nodes_scaled,'gx')
    % plot(x,y,'.k')
    % plot(x(sig_t>0),y(sig_t>0),'or')
end
function [m_dist] = get_mahalanobis_distance(nodes,tris,sig_new,sig_t,G_post)
% p-value: d^2 follows the chi-squared distribution with n degrees of freedom, 
% where n is the number of dimensions of the normal distribution. n = ?

% UNWEIGHTED
% m_dist = ((sig_t-sig_new)'/(G_post)*(sig_t-sig_new))).^0.5;

% % WEIGHTED - SIMPLE   
% x = nodes(:,1); y = nodes(:,2);
% md_wts = get_node_weights(x,y,tris);
% wtd_diff = md_wts.*(sig_t-sig_new);
% m_dist = (wtd_diff'/G_post*wtd_diff)^0.5;

% WEIGHTED - MASS MATRIX
% mass matrix
nn = length(nodes);
[pM,~,Mindex] = Premass2D(nodes,tris,0);
mm = sparse(Mindex(:,1),Mindex(:,2),pM*ones(nn,1),nn,nn);
sum(mm(:))
mm = mm/(sum(mm(:))/nn); % div by mean area per node
% md
% diff = (sig_t-sig_new);
% m_dist = (diff'*(mm\G_post)*diff)^0.5;
m_dist = ((sig_t-sig_new)'/(G_post*mm)*(sig_t-sig_new)).^0.5;
% m_dist
% m_dist/nn
m_dist = m_dist/nn;

   

end
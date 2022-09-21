function [re] = get_relative_error(nodes,tris,sig,sig_t)

% % UNWEIGHTED
% % rel_err = norm(sig_new-sig_t)/norm(sig_t);
% r_e = norm(get_truncated_sigma(sig_new)-sig_t)/norm(sig_t);

% % WEIGHTED - SIMPLE
% x = nodes(:,1); y = nodes(:,2);
% re_wts = get_node_weights(x,y,tris);
% norm_diff = (sum(re_wts.*((sig-sig_t).^2)))^.5;
% norm_truth = (sum(re_wts.*(sig_t.^2)))^.5;
% re = norm_diff/norm_truth;

% WEIGHTED - MASS MATRIX
% mass matrix
nn = length(nodes);
[pM,~,Mindex] = Premass2D(nodes,tris,0);
mm = sparse(Mindex(:,1),Mindex(:,2),pM*ones(nn,1),nn,nn);
mm = mm/(sum(mm(:)/nn));
% re
diff = (sig-sig_t);
re = (diff'*(mm*diff))/(sig_t'*(mm*sig_t));
% re

end

function [tv] = get_total_variation(nodes,tris,G_post)

% % WEIGHTED - SIMPLE
% % node weights for total variation
% x = nodes(:,1); y = nodes(:,2);
% [tv_wts] = get_node_weights(x,y,tris);
% % sum of weighted variances from posterior covariance matrix
% % similar to trace(G_post)
% tv = sum(tv_wts.*diag(G_post));

% WEIGHTED - MASS MATRIX
% mass matrix
nn = length(nodes);
[pM,~,Mindex] = Premass2D(nodes,tris,0);
mm = sparse(Mindex(:,1),Mindex(:,2),pM*ones(nn,1),nn,nn);
mm = mm/(sum(mm(:)/nn));
% tv
tv = sum(diag(G_post*mm));
tv = tv/nn;
% tv

end

% OLD FORMULATION
% function [tv_mean,tv_stdv] = get_total_variation(nodes,tris,x,y,sig_new,sigma_samples)
%     %  first calculate mesh element areas
%     tri_areas = zeros([length(tris) 1]);
%     for i = 1:length(tri_areas)
%         x1 = nodes(tris(i,1),1); y1 = nodes(tris(i,1),2);
%         x2 = nodes(tris(i,2),1); y2 = nodes(tris(i,2),2);
%         x3 = nodes(tris(i,3),1); y3 = nodes(tris(i,3),2);
%         tri_areas(i) = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
%     end
%     %  then gradient and total variance at each element (relative to mean) 
%     %  for each posterior sample
%     tot_vars = zeros([size(sigma_samples,2) 1]);
%     for ii=1:size(sigma_samples,2)
%         [dFx,dFy] = trigradient(tris,x,y,sigma_samples(:,ii)-sig_new,'face'); % at faces!
%         tot_vars(ii) = sum((abs(dFx)+abs(dFy)).*tri_areas)/sum(tri_areas);
%     end
%     tv_mean = mean(tot_vars);
%     tv_stdv = std(tot_vars);
% end




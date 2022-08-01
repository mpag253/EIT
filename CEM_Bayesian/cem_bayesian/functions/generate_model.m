function [model_vars] = generate_model(model_params)
%generate_model [Summary of this function goes here]
%   [Detailed explanation goes here]

% UNPACK
[params_m,params_e] = deal(model_params{:});
[tgeom,Nodes,Tris,full_bdy] = deal(params_m{:});
[n_elec, w_elec] = deal(params_e{:});


% GEOMETRY & MESH

% All nodes
x = Nodes(:,1); y = Nodes(:,2);                   % cartesian coords of nodes
theta = atan2(y,x);                                  % polar angles of nodes
nn = length(Nodes);                               % # nodes

% All boundary nodes
bdy_indx = find(full_bdy);                        % indices of boundary nodes
nn_bdy = sum(full_bdy);                           % # nodes on boundary
x_bdy = x(full_bdy); y_bdy = y(full_bdy);         % boundary x&y
theta_bdy = atan2(y_bdy,x_bdy);                   % polar angles of boundary nodes
[~,ss_b] = sort(theta_bdy);                       % indices when theta-sorted
% bdy_elems = [ss_b [ss_b(2:end); ss_b(1)]];        % boundary elements
bdy_elems = freeBoundary(tgeom);
% % periodic theta
% theta_bdy_p = [theta_bdy(ss_b)-2*pi; theta_bdy(ss_b); theta_bdy(ss_b)+2*pi];  
% ss_b_p = [ss_b;ss_b;ss_b];


% ELECTRODES - ANGLE BASED
% 
% % Original Electrodes
% theta_elec = linspace(-0*pi,2*pi,n_meas_pts+1)';         % polar theta array
% theta_elec = theta_meas(1:end-1);                    % for meas points
% elec_pts = .99*[interp1(theta_bdy_p,x_bdy(ss_b_p),theta_meas-pi), ...
%                 interp1(theta_bdy_p,y_bdy(ss_b_p),theta_meas-pi)];
% 
% % Define nodes in each electrode
% in_elec=zeros(nn_bdy,n_elec);
% % n_per_elec=zeros(n_elec,1);
% for ii=1:n_elec
%     theta_min = theta_elec(ii)-s_elec/2-eps;
%     theta_max = theta_elec(ii)+s_elec/2+eps;
%     if theta_min < -pi
%         theta_min = theta_min + 2*pi;
%         in_elec(:,ii) = theta_bdy >= theta_min | theta_bdy <= theta_max; 
%     elseif theta_max > pi
%         theta_max = theta_max - 2*pi;
%         in_elec(:,ii) = theta_bdy >= theta_min | theta_bdy <= theta_max; 
%     else
%         in_elec(:,ii) = theta_bdy >= theta_min & theta_bdy <= theta_max; 
%     end
% %     n_per_elec(ii)=sum(in_elec(:,ii));
% end
% in_elec=logical(in_elec);
% n_per_elec = sum(in_elec,1)';


% ELECTRODES - PERIMETER BASED

% Equidistant Electrodes
% identify boundary coords
bdy_coords = [x_bdy(ss_b) y_bdy(ss_b)];
% cumulative perimeter for each boundary node
bdy_lens = sum((bdy_coords - bdy_coords([end,1:end-1],:)).^2,2).^.5;
bdy_cper = cumsum(bdy_lens);
bdy_len = bdy_cper(end);
% find starting point (sternum...)
e0_theta = pi/2; % important: between -pi and pi
e0_cper = interp1(theta_bdy(ss_b),bdy_cper,e0_theta,'spline');
% equidistant array
bdy_elec = e0_cper + ((bdy_len/n_elec)*[0:n_elec-1])';
bdy_elec(bdy_elec>bdy_len) = bdy_elec(bdy_elec>bdy_len) - bdy_len;
bdy_elec(bdy_elec<0) = bdy_elec(bdy_elec>bdy_len) + bdy_len;
% interpolate for electrode locations
elec_scale = 1.00;
elec_pts = elec_scale*[interp1(bdy_cper,x_bdy(ss_b),bdy_elec,'spline'), ...
                       interp1(bdy_cper,y_bdy(ss_b),bdy_elec,'spline')];
% theta_elec = atan2(elec_pts(:,2),elec_pts(:,1));
% theta_elec(theta_elec<0) = theta_elec(theta_elec<0) + 2*pi;

% Define nodes in each electrode
in_elec=zeros(nn_bdy,n_elec);
% array to unsort boundary
un_ss_b(ss_b) = 1:length(ss_b);
% get boundary nodes in each electrode
for ii=1:n_elec
    e_min = bdy_elec(ii)-w_elec/2-eps;
    e_max = bdy_elec(ii)+w_elec/2+eps;
    if e_min < 0
        e_min = e_min + bdy_len;
        in_elec(:,ii) = (bdy_cper(un_ss_b)>=e_min) | (bdy_cper(un_ss_b)<=e_max); 
    elseif e_max > bdy_len
        e_max = e_max - bdy_len;
        in_elec(:,ii) = (bdy_cper(un_ss_b)>=e_min) | (bdy_cper(un_ss_b)<=e_max); 
    else
        in_elec(:,ii) = (bdy_cper(un_ss_b)>=e_min) & (bdy_cper(un_ss_b)<=e_max); 
    end
end
in_elec=logical(in_elec);
% number of nodes per electrode
n_per_elec = sum(in_elec,1)';





% Pack all variables
vars_all = {x,y,theta,nn};
vars_bdy = {x_bdy,y_bdy,theta_bdy,nn_bdy,bdy_indx,ss_b,bdy_elems}; %,theta_bdy_p};
vars_elec = {elec_pts, in_elec, n_per_elec}; %meas_op};
model_vars = {vars_all,vars_bdy,vars_elec};


end

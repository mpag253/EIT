function [model_vars] = generate_model(geometry,input_vars)
%generate_model [Summary of this function goes here]
%   [Detailed explanation goes here]

% UNPACK
[tgeom,Nodes,Tris,full_bdy] = deal(geometry{:});
[def_rob,n_meas_pts] = deal(input_vars{:});

% GEOMETRY & MESH

% All nodes
x = Nodes(:,1); y = Nodes(:,2);                   % cartesian coords of nodes
theta = atan2(y,x);                                  % polar angles of nodes
nn = length(Nodes);                               % # nodes

% All boundary nodes
bdy_indx = find(full_bdy);                        % indices of boundary nodes
nn_bdy = sum(full_bdy);                             % # nodes on boundary
x_bdy = x(full_bdy); y_bdy = y(full_bdy);
theta_bdy = atan2(y_bdy,x_bdy);       % polar angles of boundary nodes
[~,ss_b] = sort(theta_bdy);                       % indices when theta-sorted
% periodic theta neumann
theta_bdy_p = [theta_bdy(ss_b)-2*pi; theta_bdy(ss_b); theta_bdy(ss_b)+2*pi];  
ss_b_p = [ss_b;ss_b;ss_b];

% Robin boundary condition
% make this function of theta
% bdy_rob = full_bdy & abs(x)<=xtol_rob & y<0;      % selection of robin BC
bdy_rob = zeros(length(x),1);
for k=1:length(def_rob)
    for j=1:3  % to capture periodic bounds
        bdy_rob = bdy_rob + (full_bdy & theta>=def_rob{k}(1)+(j-2)*2*pi & ...
                                        theta<=def_rob{k}(2)+(j-2)*2*pi);
    end
end
bdy_rob = logical(bdy_rob);
bdy_rob_indx = find(bdy_rob);                     % indices of robin BC
nn_rob = sum(bdy_rob);                            % # nodes with robin BC
[~,ss_rob] = sort(atan2(y(bdy_rob),x(bdy_rob)));  % indices when theta-sorted

% Neumann boundary condition
bdy_neu = full_bdy & ~bdy_rob;                    % neumann BC if not robin
bdy_neu_indx = find(bdy_neu);                     % indices of neumann BC
nn_neu = sum(bdy_neu);                            % # nodes with neumann BC
[~,ss_neu] = sort(atan2(y(bdy_neu),x(bdy_neu)));  % indices when theta-sorted
theta_neu_c = atan2(y(bdy_neu),x(bdy_neu));       % theta of neumann nodes    (_c=circle?)
theta_neu_p = [theta_neu_c-2*pi; ...
               theta_neu_c; ...
               theta_neu_c+2*pi];                 % periodic theta neumann

% ELECTRODES
% % Original
% theta_meas = linspace(-0*pi,2*pi,n_meas_pts+1)';         % polar theta array
% theta_meas = theta_meas(1:end-1);                    % for meas points
% meas_pts = .99*[interp1(theta_bdy_p,x_bdy(ss_b_p),theta_meas-pi), ...
%                 interp1(theta_bdy_p,y_bdy(ss_b_p),theta_meas-pi)];
% Equidistant
% identify boundary nodes
bdy_nods = [x_bdy(ss_b) y_bdy(ss_b)];  
% cumulative perimeter for each boundary node
bdy_lens = sum((bdy_nods - bdy_nods([end,1:end-1],:)).^2,2).^.5;
bdy_cper = cumsum(bdy_lens);
bdy_perim = bdy_cper(end);
%find starting point (sternum...)
init_pt_theta = -pi;
bdy_cper_p = [bdy_cper-bdy_perim;bdy_cper;bdy_cper+bdy_perim];
init_pt_cper = interp1(theta_bdy_p,bdy_cper_p,init_pt_theta);
% equidistant array
bdy_meas = init_pt_cper + (bdy_cper(end)/(n_meas_pts)*[0:n_meas_pts-1])';
bdy_meas(bdy_meas>bdy_perim) = bdy_meas(bdy_meas>bdy_perim) - bdy_perim;
% interpolate for electrode locations
meas_pts = .99*[interp1(bdy_cper_p,x_bdy(ss_b_p),bdy_meas), ...
                interp1(bdy_cper_p,y_bdy(ss_b_p),bdy_meas)];
theta_meas = atan2(meas_pts(:,2),meas_pts(:,1));
theta_meas(theta_meas<0) = theta_meas(theta_meas<0) + 2*pi;



%  Pack geometry and mesh veriables
vars_all = {x,y,theta,nn};
vars_bdy = {x_bdy,y_bdy,theta_bdy,bdy_indx,nn_bdy,ss_b,theta_bdy_p};  % bdy_indx
vars_rob = {bdy_rob,bdy_rob_indx,nn_rob,ss_rob};
vars_neu = {bdy_neu,bdy_neu_indx,nn_neu,ss_neu,theta_neu_c,theta_neu_p};
vars_mea = {theta_meas,meas_pts};


%%%%%%%%%%%%%%%%%%%%%%%%% Meas & Bdy Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neumann boundary condition   
neu_op = zeros(nn_neu,nn);                      % operator for neumann nodes
neu_op(1:nn_neu,bdy_neu_indx) = eye(nn_neu);    % 1s identify neumann nodes

% Robin boundary condition
rob_op = zeros(nn_rob,nn);                      % operator for robin nodes
rob_op(1:nn_rob,bdy_rob_indx) = eye(nn_rob);    % 1s identify robin nodes

% Measurement electrodes
meas_op = zeros(n_meas_pts,nn);               % operator for meas nodes
for ii = 1:n_meas_pts 
    % find the triangle enclosing the measurement point    
    [ti,barycoords] = pointLocation(tgeom,meas_pts(ii,:));
    % store the nodal weights (barycentric) in the operator
    meas_op(ii,Tris(ti,:)) = barycoords;
end

% Pack operators
vars_op = {neu_op, rob_op, meas_op};

% Pack all variables
model_vars = {vars_all,vars_bdy,vars_rob,vars_neu,vars_mea,vars_op};


end

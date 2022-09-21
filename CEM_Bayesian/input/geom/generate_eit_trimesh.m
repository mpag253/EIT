clear, clc, close all
addpath('../../cem_bayesian/distmesh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  INPUTS

pca_id = 'pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A';
version = 'truth';
% version = 'truth_inv';
% version = 'sample_mean';
% version = 'predicted';

% Version
if pca_id(end) == 'A', subject_id = 'H5977';
elseif pca_id(end) == 'B', subject_id = 'AGING043';
elseif pca_id(end) == 'C', subject_id = 'H7395';
elseif pca_id(end) == 'D', subject_id = 'AGING014';
elseif pca_id(end) == 'E', subject_id = 'AGING053';
end


if strcmp(version,'truth')
    version_id = strcat('truth_', subject_id);
    mesh_size_t = 1.5; %2.5;  % torso
    mesh_size_i = 6; %5;  % interior
    refine_grad = 0.3;  % slow=0.1, fast=0.3
    path = strcat(pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat('mesh_seeds_torso_',sprintf('%04d',10*mesh_size_t),'.csv');
    pv_file_l = strcat('mesh_seeds_lung_left_',sprintf('%04d',10*mesh_size_i),'.csv');
    pv_file_r = strcat('mesh_seeds_lung_right_',sprintf('%04d',10*mesh_size_i),'.csv');
    mesh_lungs = true;

elseif strcmp(version,'truth_inv')
    version_id = strcat('truth_', subject_id);
    mesh_size_t = 2;  % torso
    mesh_size_i = 8;  % interior
    refine_grad = 0.3;  % slow=0.1, fast=0.3
    path = strcat(pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat('mesh_seeds_torso_',sprintf('%04d',10*mesh_size_t),'.csv');
    mesh_lungs = false;

elseif strcmp(version,'sample_mean')
    version_id = 'sample_mean';
    mesh_size_t = 2;  % torso
    mesh_size_i = 8;  % interior
    refine_grad = 0.3;  % slow=0.1, fast=0.3
    path = strcat(pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat('mesh_seeds_torso_',sprintf('%04d',10*mesh_size_t),'.csv');
    mesh_lungs = false;

elseif strcmp(version,'predicted')
    version_id = strcat('predicted_', subject_id);
    mesh_size_t = 2;  % torso
    mesh_size_i = 8;  % interior
    refine_grad = 0.3;  % slow=0.1, fast=0.3
    path = strcat(pca_id,'/',version_id,'/mesh_seeds/');
    pv_file_t = strcat('mesh_seeds_torso_',sprintf('%04d',10*mesh_size_t),'.csv');
    mesh_lungs = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% READ DATA

pv_t = readmatrix(strcat(path,pv_file_t));
if mesh_lungs == true
    pv_l = readmatrix(strcat(path,pv_file_l));
    pv_r = readmatrix(strcat(path,pv_file_r));
    pv_all = [pv_t;pv_l;pv_r];
else
    pv_all = pv_t;
end

mesh_area = polyarea(pv_t(:,1),pv_t(:,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATE TRIMESH

bbox = 1.01*[min(pv_t(:,1)) min(pv_t(:,2)); max(pv_t(:,1)) max(pv_t(:,2))];
% [nodes,tris]=distmesh2d(@dpoly,@huniform,mesh_size,bbox,pv,pv);
% [nodes,tris]=distmesh2d(@dpoly,@huniform,mesh_size_i,bbox,pv_all,pv_t);
fd = @(p) dpoly(p, pv_t); %@dpoly;
fh = @(p) min(mesh_size_t-refine_grad*dpoly(p, pv_t), mesh_size_i);
h0 = min(mesh_size_t,mesh_size_i);
fix = pv_all;
[nodes,tris]=distmesh2d(fd,fh,h0,bbox,fix);

tgeom = triangulation(tris,nodes);
% tgeom = delaunayTriangulation(Nodes);
% Tris = tgeom.ConnectivityList;

% boundary_indices = boundary(nodes(:,1), nodes(:,2), shrink_factor);
% boundary_operator(boundary_indices) = true;
% boundary_operator = boundary_operator';

mesh_specifics = strcat(sprintf('%04d',10*mesh_size_i),'_',sprintf('%04d',10*mesh_size_t));
save_path = strrep(path, 'mesh_seeds/', strcat('trimesh/trimesh_',mesh_specifics,'/'));
if ~exist(save_path, 'dir'), mkdir(save_path), end
save(strcat(save_path, 'trimesh.mat'),'tgeom')
save_file_nods = strcat('trimesh_nodes.csv');
save_file_tris = strcat('trimesh_nodes.csv');
writematrix(tgeom.Points, strcat(save_path,save_file_nods))
writematrix(tgeom.ConnectivityList, strcat(save_path,save_file_tris))

disp("Done.")

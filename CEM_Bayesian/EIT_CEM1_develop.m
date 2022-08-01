
clear, clc, close all
addpath('cem_bayesian')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESH
pca_id = 'pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A';
fwd_mesh = 'truth_H5977/trimesh/trimesh_0100_0050';
% inv_mesh = 'sample_mean/trimesh/trimesh_0150_0050';
inv_mesh = 'truth_H5977/trimesh/trimesh_0100_0050';
% inv_mesh = 'truth_H5977/trimesh/trimesh_0150_0050_inv';

% ELECTRODES & PATTERNS
n_elec = 16;                    % # electrodes
n_patt = n_elec;                % # measurement patterns
z_elec = 0.001;                 % contact impedance
w_elec = 20;                    % width of electrode

% DATA
noise = .1 %1.;                        % noise level [%]
sig_0 = 1;                    % (...)

% PRIORS
n_priors = 9;                   % # priors to plot
l_prcorr = 8e1;                 % prior correlation length
prior_type = 'p0';              % p0=gauss, p1=mean+gauss, p2=mean+gauss*cov


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_m = {pca_id,fwd_mesh,inv_mesh};
params_e = {n_elec,n_patt,z_elec,w_elec};
params_d = {noise,sig_0};
params_p = {prior_type,n_priors,l_prcorr};
eit_params = {params_m,params_e,params_d,params_p};
run_eit_cem_bayes(eit_params)


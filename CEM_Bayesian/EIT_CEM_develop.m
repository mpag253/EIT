
clear, clc, close all
addpath('cem_bayesian')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESH
pca_id = 'pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-A';
[fwd_mesh,fwd_s1,fwd_s2] = deal('truth_H5977','0060','0015');
% [inv_mesh,inv_s1,inv_s2] = deal('predicted_H5977','0080','0020');  
[inv_mesh,inv_s1,inv_s2] = deal('sample_mean','0080','0020');
% [inv_mesh,inv_s1,inv_s2] = deal('truth_H5977','0080','0020');
% [inv_mesh,inv_s1,inv_s2] = deal('truth_H5977','0060','0015');

% pca_id = 'pca_LRT_S_Mfull_N80_R-AGING001-EIsupine_LOO-B';
% [fwd_mesh,fwd_s1,fwd_s2] = deal('truth_AGING043','0060','0015');
% [inv_mesh,inv_s1,inv_s2] = deal('predicted_AGING043','0080','0020'); 

% ELECTRODES & PATTERNS
n_elec = 16;                    % # electrodes
z_elec = 0.001;                 % contact impedance
w_elec = 20;                    % width of electrode

% DATA
noise = 1.0;                    % noise level [%]
conditions = [];%["RA",];                % diseased conditions ("RA","RP","LA","LP")

% PRIORS
prior_type = 'p0';              % p0=gauss, p1=mean+gauss, p2=mean+gauss*cov
l_prcorr = 8e1;                 % prior correlation length

% RUN
make_figs = [103,105,106];  % [101,102,103,104,105,106]; 
save_figs = [];
save_data = false;
batch_name = 'develop';
num = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_m = {pca_id,fwd_mesh,fwd_s1,fwd_s2,inv_mesh,inv_s1,inv_s2};
params_e = {n_elec,z_elec,w_elec};
params_d = {noise,conditions};
params_p = {prior_type,l_prcorr,};
eit_params = {params_m,params_e,params_d,params_p};
run_params = {make_figs,save_figs,save_data,batch_name,num2str(num)};
run_eit_cem_bayes(eit_params,run_params)

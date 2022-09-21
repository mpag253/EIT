
clear, clc, close all
addpath('cem_bayesian')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Batch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SETUP
batch_name = 'batch_manuscript_cases';
% batch_name = 'batch_temporal_demo_2';
make_figs = [103,106];  % [101,102,103,104,105,106]; 
save_figs = [103,106];
save_data = true;
fname = ['output/metrics/metrics_',batch_name,'.xlsx'];
if isfile(fname), delete(fname); end

% RUN
batch_file = ['input/batches/input_',batch_name,'.xlsx'];
inputs = readcell(batch_file,'NumHeaderLines',1);
for i = 1:size(inputs,1)
    do = inputs(i,1);
    if do{1}

        % Allocate parameters
        inp = inputs(i,2:end);
        num = inp{1};
        params_m = inp(2:8);    % {pca_id,fwd_mesh,fwd_s1,fwd_s2,inv_mesh,inv_s1,inv_s2}
        params_e = inp(9:11);   % {n_elec,z_elec,w_elec};
        params_d = inp(12:14);  % {noise,conditions};
        params_p = inp(15:16);  % {prior_type,l_prcorr};
        params_d{2} = string(split(params_d{2},','));  % format conditions
        eit_params = {params_m,params_e,params_d,params_p};
        run_params = {make_figs,save_figs,save_data,batch_name,num2str(num)};
    
        % Run
        run_eit_cem_bayes(eit_params,run_params)
        close("all")
        
    end
end




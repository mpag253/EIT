function [section_data] = get_sections(x,y,x_s,y_s,sig_t,sig_new,S_post,sig_samps)

    % Truth
    F = scatteredInterpolant(x,y,sig_t);
    z_t = F(x_s,y_s);
    
    % MAP estimate
    F = scatteredInterpolant(x,y,sig_new);         
    z_s = F(x_s,y_s);
    
    % Posterior standard devaition
    F = scatteredInterpolant(x,y,S_post);          
    z_s_sd = F(x_s,y_s);
    
    % Posterior samples
    n_samps = size(sig_samps,2);
    z_s_samps = cell(n_samps);
    for i = 1:n_samps
        F = scatteredInterpolant(x,y,sig_samps(:,1));
        z_s_samps{i} = F(x_s,y_s);
    end
    
    % Confidence interval
    z_s_ci = {z_s+2.576*z_s_sd, z_s-2.576*z_s_sd};

    % Pack output
    section_data = {z_t,z_s,z_s_samps,z_s_ci};

end
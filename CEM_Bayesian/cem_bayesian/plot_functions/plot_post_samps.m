function [] = plot_post_samps(pt,fg,sp,tris,x,y,sig_samps,clims,savefile)

    for ii=1:size(sig_samps,2)

        % Plot
        % sig_samp_LS = get_level_set(sig_samps(:,ii));
        plot_sigma(pt,fg,sp,tris,x,y,sig_samps(:,ii),clims)
        % plot_sigma('Posterior Samples',1,[2,2,4],tris,x,y,sig_samp_LS,clims1)
          
        if savefile
            % Capture the plot as an image 
            frame = getframe(figure(1)); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256);
    
            % Write to the GIF File 
            if ii == 1, imwrite(imind,cm,savefile,'gif', 'Loopcount',inf); 
            else, imwrite(imind,cm,savefile,'gif','WriteMode','append'); end 
        end

    end
end
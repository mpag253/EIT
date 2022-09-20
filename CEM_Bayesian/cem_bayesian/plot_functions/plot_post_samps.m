function [] = plot_post_samps(pt,fg,sp,tris,x,y,sig_samps,clims,save_dir)

    for ii=1:size(sig_samps,2)

        % Plot
        % sig_samp_LS = get_level_set(sig_samps(:,ii));
        plot_tau(pt,fg,sp,tris,x,y,sig_samps(:,ii),clims)
        % plot_sigma('Posterior Samples',1,[2,2,4],tris,x,y,sig_samp_LS,clims1)
        set(gcf,'Color',[1,1,1]); set(gca,'Color',[1,1,1]);
          
        if save_dir
            % Capture the plot as an image 
            frame = getframe(figure(fg)); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256);

            % Write first 10 frames to images
            if ii < 11
                save_file = [save_dir,'figure_106_',num2str(ii,'%02.f'),'.png'];
                imwrite(im,save_file)
            end
    
            % Write to the GIF File
            save_file = [save_dir,'figure_106_animated.gif'];
            if ii == 1, imwrite(imind,cm,save_file,'gif', 'Loopcount',inf); 
            else, imwrite(imind,cm,save_file,'gif','WriteMode','append'); end 
        end

    end
end
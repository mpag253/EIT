
root = 'C:\Users\mipag_000\OneDrive - The University of Auckland\Documents\EIT_Project\EIT\CEM_Bayesian\';
path = 'output\figures\batch_temporal_demo_2\';

for i = 1:85

    fname = [root path 'batch_temporal_demo_2-' num2str(i) '\figure_103.png'];
    im = imread(fname);
    [imind,cm] = rgb2ind(im,256);

    save_file = [root path 'batch_temporal_demo_2.gif'];
    delay = 0.02;
    if i == 1, imwrite(imind,cm,save_file,'gif','DelayTime',delay,'Loopcount',inf); 
    else, imwrite(imind,cm,save_file,'gif','DelayTime',delay,'WriteMode','append'); end 
    
end

for i = 85:-1:1

    fname = [root path 'batch_temporal_demo_2-' num2str(i) '\figure_103.png'];
    im = imread(fname);
    [imind,cm] = rgb2ind(im,256);

    save_file = [root path 'batch_temporal_demo_2.gif'];
    delay = 0.02;
    imwrite(imind,cm,save_file,'gif','DelayTime',delay,'WriteMode','append');
    
end








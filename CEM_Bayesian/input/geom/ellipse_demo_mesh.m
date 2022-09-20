addpath('distmesh')

% Example: (Ellipse)
fd=@(p) p(:,1).^2/4^2+p(:,2).^2/3^2-1;
[p,t]=distmesh2d(fd,@huniform,0.1,[-2,-1;2,1],[]);
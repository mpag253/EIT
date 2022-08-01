function phi=gaussian(centre,width,height,x)
% phi=gaussian(centre,width,height,x)
%
if min(size(x))==1
    f=exp(-((x-centre).^2/(width^2)));
    phi=height*f/max(f);
elseif min(size(x))==2
    if length(width)==1 
        width=[width,width]; 
    end

    f=exp(-(((x(:,1)-centre(1)).^2)/(width(1)^2)+...
        ((x(:,2)-centre(2)).^2)/(width(2)^2)));
    phi=height*f/max(f);
else
    if length(width)==1
        width=[width,width,width]; 
    end
    f=exp(-(((x(:,1)-centre(1)).^2)/(width(1)^2)+...
        ((x(:,2)-centre(2)).^2)/(width(2)^2)+...
        ((x(:,3)-centre(3)).^2)/(width(3)^2)));
    phi=height*f/max(f);
end
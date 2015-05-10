function [ img ] = AddImages( imgs )
%[ img ] = AddImages( imgs )


[x,y,z] = size(imgs);
img = zeros(x,y);

for n = 1:z
    img = img + imgs(:,:,n);
end

end


function [ azi_avg ] = azimuthal_average( A, center_point )
%UNTITLED Summary of this function goes here
%  Compute the average of a matrix A azimuthally


[m,n] = size(A);

if nargin < 2
    center_point = [round(m/2),round(n/2)];
end


[i,j] = meshgrid(1:m,1:n);
i = i(:); j = j(:);

dist = round(sqrt((i-center_point(1)).^2 + (j-center_point(2)).^2));
[dist,y] = sort(dist);


hh = hist(dist,max(dist)+1);

vec = A(:);
vec = vec(y); % sort the same way as dist

ini = 2;

azi_avg(1:max(dist))=0;

for k = 1:max(dist)

   index = [ini:ini+hh(k+1)-1]; 
   RoI = vec(index);
   azi_avg(k) = mean(RoI(abs(RoI)>0));
   ini = max(index)+1;

end

azi_avg(isnan(azi_avg)) = 0;












end


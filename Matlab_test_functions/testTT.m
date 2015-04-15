clear all
clc
close all

load('dOTF_phase.mat')

lambda = 0.6;
k = (2*pi)/lambda;
mask_radius = 5;

imagesc(dOTF_phase)

p2 = pickPoint;
p1 = pickPoint;
p3 = pickPoint;
p7 = pickPoint;
p8 = pickPoint;
p9 = pickPoint;

points = cell(1,3);
points{1} = p2;
points{2} = p1;
points{3} = p3;
points{4} = p7;
points{5} = p8;
points{6} = p9;
                        
[sizey,sizex] = size(dOTF_phase);
xx = linspace(1,sizex,sizex);
yy = linspace(1,sizey,sizey);
[X,Y] = meshgrid(xx,yy);

R2 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
M2 = R2<=mask_radius;
[tip2,tilt2] = calctiptiltdOTF(M2.*dOTF_phase*k,p2,[3,3]);


R1 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
M1 = R1<=mask_radius;
[tip1,tilt1] = calctiptiltdOTF(M1.*dOTF_phase*k,p1,[3,3]);


R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
M3 = R3<=mask_radius;
[tip3,tilt3] = calctiptiltdOTF(M3.*dOTF_phase*k,p3,[3,3]);


R7 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
M7 = R7<=mask_radius;
[tip7,tilt7] = calctiptiltdOTF(M7.*dOTF_phase*k,p7,[3,3]);

R8 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
M8 = R8<=mask_radius;
[tip8,tilt8] = calctiptiltdOTF(M8.*dOTF_phase*k,p8,[3,3]);

R9 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
M9 = R9<=mask_radius;
[tip9,tilt9] = calctiptiltdOTF(M9.*dOTF_phase*k,p9,[3,3]);

PTTpos = zeros(37,3);
PTTpos(1,2) = tip1;
PTTpos(1,3) = tilt1;
PTTpos(2,2) = tip2;
PTTpos(2,3) = tilt2;
PTTpos(3,2) = tip3;
PTTpos(3,3) = tilt3;
PTTpos(7,2) = tip7;
PTTpos(7,3) = tilt7;
PTTpos(8,2) = tip8;
PTTpos(8,3) = tilt8;
PTTpos(9,2) = tip9;
PTTpos(9,3) = tilt9;

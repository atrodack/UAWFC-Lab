clear all
clc
close all

load('dOTF_phase.mat')

lambda = 0.6;
k = (2*pi)/lambda;
mask_radius = 5;

[sy,sx] = size(dOTF_phase);


[fx,fy] = gradient(dOTF_phase .* k,10,10);
imagesc(fx)


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
                        
[sizey,sizex] = size(fx);
xx = linspace(1,sizex,sizex);
yy = linspace(1,sizey,sizey);
[X,Y] = meshgrid(xx,yy);



R2 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
M2 = R2<=mask_radius;
seg2tilt = M2.*fx;
seg2tip = M2.*fy;
seg2Proptilt = seg2tilt(abs(seg2tilt)>0);
seg2Proptip = seg2tip(abs(seg2tip)>0);
tilt2 = mean(seg2Proptilt);
tip2 = mean(seg2Proptip);


R1 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
M1 = R1<=mask_radius;
seg1tilt = M1.*fx;
seg1tip = M1.*fy;
seg1Proptilt = seg1tilt(abs(seg1tilt)>0);
seg1Proptip = seg1tip(abs(seg1tip)>0);
tilt1 = mean(seg1Proptilt);
tip1 = mean(seg1Proptip);

R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
M3 = R3<=mask_radius;
seg3tilt = M3.*fx;
seg3tip = M3.*fy;
seg3Proptilt = seg3tilt(abs(seg3tilt)>0);
seg3Proptip = seg3tip(abs(seg3tip)>0);
tilt3 = mean(seg3Proptilt);
tip3 = mean(seg3Proptip);

R7 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
M7 = R7<=mask_radius;
seg7tilt = M7.*fx;
seg7tip = M7.*fy;
seg7Proptilt = seg7tilt(abs(seg7tilt)>0);
seg7Proptip = seg7tip(abs(seg7tip)>0);
tilt7 = mean(seg7Proptilt);
tip7 = mean(seg7Proptip);

R8 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
M8 = R8<=mask_radius;
seg8tilt = M8.*fx;
seg8tip = M8.*fy;
seg8Proptilt = seg8tilt(abs(seg8tilt)>0);
seg8Proptip = seg8tip(abs(seg8tip)>0);
tilt8 = mean(seg8Proptilt);
tip8 = mean(seg8Proptip);

R9 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
M9 = R9<=mask_radius;
seg9tilt = M9.*fx;
seg9tip = M9.*fy;
seg9Proptilt = seg9tilt(abs(seg9tilt)>0);
seg9Proptip = seg9tip(abs(seg9tip)>0);
tilt9 = mean(seg9Proptilt);
tip9 = mean(seg9Proptip);

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
PTTpos = PTTpos .* 1e3;

subplot(1,2,1);
imagesc(seg1tilt+seg2tilt+seg3tilt+seg7tilt+seg8tilt+seg9tilt);
subplot(1,2,2)
imagesc(seg1tip+seg2tip+seg3tip+seg7tip+seg8tip+seg9tip);







% phasegrabbingmask

clear all;
clc;
close all;

cd('/home/lab/Desktop/dOTF_saved_structures')
load('dOTFstructure_data_run_20152101.mat');
cd /home/lab/Desktop

mask_radius = 5; %pixels

dotf = dOTF.storeddOTF{1};
dOTF.Phase = angle(dotf);
% dOTF.unwrapphase('gold');
phase = dOTF.Phase;
phase = phase ./ ((2*pi)/0.6e-6);
imagesc(phase);
colormap(jet);
daspect([1,1,1]);

input('Press Enter\n');
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

[sizey,sizex] = size(phase);
xx = linspace(1,sizex,sizex);
yy = linspace(1,sizey,sizey);
[X,Y] = meshgrid(xx,yy);

R2 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
M2 = R2<=mask_radius;
seg2 = M2.*phase;
seg2piston = mean(seg2(:));

R1 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
M1 = R1<=mask_radius;
seg1 = M1.*phase;
seg1piston = mean(seg1(:));

R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
M3 = R3<=mask_radius;
seg3 = M3.*phase;
seg3piston = mean(seg3(:));

R7 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
M7 = R7<=mask_radius;
seg7 = M7.*phase;
seg7piston = mean(seg7(:));

R8 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
M8 = R8<=mask_radius;
seg8 = M8.*phase;
seg8piston = mean(seg8(:));

R9 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
M9 = R9<=mask_radius;
seg9 = M9.*phase;
seg9piston = mean(seg9(:));

clf;
imagesc(seg2+seg1+seg3+seg7+seg8+seg9);
axis xy;
daspect([1,1,1]);

PTTpos = zeros(37,3);
PTTpos(2,1) = seg2piston;
PTTpos(1,1) = seg1piston;
PTTpos(3,1) = seg3piston;
PTTpos(7,1) = seg7piston;
PTTpos(8,1) = seg8piston;
PTTpos(9,1) = seg9piston;
%Convert to base unit of micron
PTTpos = PTTpos .* 1e6;
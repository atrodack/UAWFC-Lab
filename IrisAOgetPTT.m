function [ PTT_mirror, PTTpos_mirror ] = IrisAOgetPTT( dOTF, pixelshift, lambda, overlapsegs, pokeseg, calibration_filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(calibration_filename);

k = (2*pi) / lambda;
phase = angle(dOTF);
unwrapped_phase = uwrap(phase,'unwt');
OPL = unwrapped_phase / k;


g = dOTF;
tipleft = conj(circshift(g,[0,pixelshift(2)]));
tipright = (circshift(g,[0,-pixelshift(2)]));
tiltup = (circshift(g,[-pixelshift(1),0]));
tiltdown = conj(circshift(g,[pixelshift(1),0]));

Tipdiff= tipleft .* tipright;
Tiltdiff = tiltup .* tiltdown;


[sizey,sizex] = size(phase);
xx = linspace(1,sizex,sizex);
yy = linspace(1,sizey,sizey);
[XX,YY] = meshgrid(xx,yy);

SegMasks = cell(37,1);
for n = 1:37
    if n ~= pokeseg
        R = sqrt((XX - pixel_seg_map{n}(2)).^2 + (YY - pixel_seg_map{n}(1)).^2);
        SegMasks{n,1} = double(R <= Areal_Averaging_radius);
    end
end


pistonlist = zeros(37,1);
% for n = 1:37
%     if n ~= pokeseg
%         segpiston = OPL .* SegMasks{n};
%         pistonlist(n,1) = mean(mean(segpiston));
%     end
% end
for n = 1:37
    if n ~= pokeseg
        pistonlist(n,1) = OPL(pixel_seg_map{n}(1),pixel_seg_map{n}(2));
    end
end


tiplist = zeros(37,1);
for n = 1:37
    if n ~= pokeseg
        segtilt = angle(Tipdiff .* SegMasks{n});
        tiplist(n,1) = mean(mean(segtilt));
    end
end

tiltlist = zeros(37,1);
for n = 1:37
    if n ~= pokeseg
        segtip = angle(Tiltdiff .* SegMasks{n});
        tiltlist(n,1) = mean(mean(segtip));
    end
end

piston_float = pistonlist(1);
if piston_float >= 0
    pistonlist = pistonlist - piston_float;
else
    pistonlist = pistonlist + abs(piston_float);
end

pistonlist = pistonlist/2;
% tiplist = tiplist;
% tiltlist = tiltlist;


PTTpos_mirror = horzcat(pistonlist, tiplist,tiltlist);
PTTpos_mirror(overlapsegs,:) = 0; %account for overlap region
PTT_mirror = mapSegments(PTTpos_mirror);
end


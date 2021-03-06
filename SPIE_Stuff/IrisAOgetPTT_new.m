function [ PTT_mirror, PTTpos_mirror ] = IrisAOgetPTT_new( dOTF, pixelshift, lambda, overlapsegs, pokeseg, calibration_filename, numsegs)
%[ PTT_mirror, PTTpos_mirror ] = IrisAOgetPTT( dOTF, pixelshift, lambda, overlapsegs, pokeseg, calibration_filename )
%   Detailed explanation goes here

if nargin < 7
    numsegs = 37;
    numrings = 3;
    NECO = false;
end

if numsegs == 19
    numrings = 2;
    NECO = true;
end

load(calibration_filename);

k = (2*pi) / lambda;
dx = 6.0600e-6;
phase = angle(dOTF);
unwrapped_phase = uwrap(phase,'gold');
OPL = unwrapped_phase / k;


%% Compute Tip/Tilt
g = dOTF;
tipleft = conj(circshift(g,[0,pixelshift(2)]));
tipright = (circshift(g,[0,-pixelshift(2)]));
tiltup = (circshift(g,[-pixelshift(1),0]));
tiltdown = conj(circshift(g,[pixelshift(1),0]));

Tipdiff= tipleft .* tipright;
Tipdiff = Tipdiff / (2*pixelshift(2));
Tipdiff = Tipdiff / k;
Tipdiff = Tipdiff / dx;
Tiltdiff = tiltup .* tiltdown;
Tiltdiff = Tiltdiff / (2*pixelshift(1));
Tiltdiff = Tiltdiff / k;
Tiltdiff = Tiltdiff / dx;


[sizey,sizex] = size(phase);
xx = linspace(1,sizex,sizex);
yy = linspace(1,sizey,sizey);
[XX,YY] = meshgrid(xx,yy);

SegMasks = cell(numsegs,1);
for n = 1:numsegs
    if n ~= pokeseg
        R = sqrt((XX - pixel_seg_map{n}(2)).^2 + (YY - pixel_seg_map{n}(1)).^2);
        SegMasks{n,1} = double(R <= Areal_Averaging_radius);
    end
end



tiplist = zeros(numsegs,1);
for n = 1:numsegs
    if n ~= pokeseg
        segtilt = (Tipdiff .* SegMasks{n});
        tiplist(n,1) = (mean(mean(angle(segtilt(abs(segtilt)>0)))));
    end
end



tiltlist = zeros(numsegs,1);
for n = 1:numsegs
    if n ~= pokeseg
        segtip = (Tiltdiff .* SegMasks{n});
        tiltlist(n,1) = (mean(mean(angle(segtip(abs(segtip)>0)))));
    end
end
tiplist = tiplist * 1e-2;
tiltlist = tiltlist * 1e-2;

%% Correct for Tip/Tilt and get piston
% pistonlist = zeros(numsegs,1);
% for n = 1:numsegs
%     if n ~= pokeseg
%         segpiston = OPL .* SegMasks{n};
%         segpiston = fftshift(circshift(segpiston,1 - [pixel_seg_map{n}(1),pixel_seg_map{n}(2)]));
%         corrector = flipud(fliplr(segpiston));
%         tt_removed = (segpiston + corrector) / 2;
%         pistonlist(n) = mean(mean(tt_removed(abs(tt_removed)>0)));
%     end
% end




pistonlist = zeros(numsegs,1);
for n = 1:numsegs
    if n ~= pokeseg
        segpiston = OPL .* SegMasks{n};
        pistonlist(n,1) = mean(mean(segpiston(abs(segpiston)>0)));
    end
end



% for n = 1:numsegs
%     if n ~= pokeseg
%         pistonlist(n,1) = OPL(pixel_seg_map{n}(1),pixel_seg_map{n}(2));
%     end
% end

maxpiston = max(pistonlist(abs(pistonlist)>0));
minpiston = min(pistonlist(abs(pistonlist)>0));
piston_float = (abs(maxpiston)+abs(minpiston))/2;
if abs(maxpiston) >= abs(minpiston)
    pistonlist(abs(pistonlist)>0) = pistonlist(abs(pistonlist)>0) - piston_float;
else
    pistonlist(abs(pistonlist)>0) = pistonlist(abs(pistonlist)>0) + piston_float;
end



% pistonlist = pistonlist * 0.65;



PTTpos_mirror = horzcat(pistonlist, tiplist,tiltlist);
PTTpos_mirror(overlapsegs,:) = 0; %account for overlap region


PTTpos_mirror(overlapsegs(1),1) = mean(PTTpos_mirror([16,6,7],1));
PTTpos_mirror(overlapsegs(3),1) = mean(PTTpos_mirror([2,7,8],1));
PTTpos_mirror(overlapsegs(2),1) = mean(PTTpos_mirror([17,19,7],1));







PTT_mirror = mapSegments(PTTpos_mirror,numrings);
end


function [ PTT_mirror, PTTpos_mirror ] = IrisAOgetPTT( dOTF, pixelshift, lambda, overlapsegs, pokeseg, calibration_filename )
%[ PTT_mirror, PTTpos_mirror ] = IrisAOgetPTT( dOTF, pixelshift, lambda, overlapsegs, pokeseg, calibration_filename )
%   Detailed explanation goes here

load(calibration_filename);

k = (2*pi) / lambda;
dx = 6.0600e-6;
phase = angle(dOTF);
unwrapped_phase = uwrap(phase,'unwt');
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

SegMasks = cell(37,1);
for n = 1:37
    if n ~= pokeseg
        R = sqrt((XX - pixel_seg_map{n}(2)).^2 + (YY - pixel_seg_map{n}(1)).^2);
        SegMasks{n,1} = double(R <= Areal_Averaging_radius);
    end
end



tiplist = zeros(37,1);
for n = 1:37
    if n ~= pokeseg
        segtilt = angle(Tipdiff .* SegMasks{n});
        tiplist(n,1) = (mean(mean(segtilt)));
    end
end
tiplist(tiplist>5e-3) = 5e-3;
tiplist(tiplist<-5e-3) = -5e-3;



tiltlist = zeros(37,1);
for n = 1:37
    if n ~= pokeseg
        segtip = angle(Tiltdiff .* SegMasks{n});
        tiltlist(n,1) = (mean(mean(segtip)));
    end
end
tiltlist(tiltlist>5e-3) = 5e-3;
tiltlist(tiltlist<-5e-3) = -5e-3;

%% Correct for Tip/Tilt and get piston
pistonlist = zeros(37,1);
% for n = 1:37
%     if n ~= pokeseg
%         segpiston = OPL .* SegMasks{n};
%         segpiston = fftshift(circshift(segpiston,1 - [pixel_seg_map{n}(1),pixel_seg_map{n}(2)]));
%         corrector = flipud(fliplr(segpiston));
%         tt_removed = (segpiston + corrector) / 2;
%         pistonlist(n) = mean(mean(tt_removed(abs(tt_removed)>0)));
%     end
% end




% pistonlist = zeros(37,1);
% for n = 1:37
%     if n ~= pokeseg
%         segpiston = OPL .* SegMasks{n};
%         pistonlist(n,1) = mean(mean(segpiston(abs(segpiston)>0)));
%     end
% end



for n = 1:37
    if n ~= pokeseg
        pistonlist(n,1) = OPL(pixel_seg_map{n}(1),pixel_seg_map{n}(2));
    end
end

maxpiston = max(pistonlist);
minpiston = min(pistonlist);
piston_float = (abs(maxpiston)+abs(minpiston))/2;
if abs(maxpiston) >= abs(minpiston)
    pistonlist = pistonlist - piston_float;
else
    pistonlist = pistonlist + piston_float;
end



pistonlist = pistonlist/2;
% tiplist = tiplist;
% tiltlist = tiltlist;


PTTpos_mirror = horzcat(pistonlist, tiplist,tiltlist);
PTTpos_mirror(overlapsegs,:) = 0; %account for overlap region
PTT_mirror = mapSegments(PTTpos_mirror);
end


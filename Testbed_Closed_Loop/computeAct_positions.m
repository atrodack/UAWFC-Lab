function [ mirror_shape ] = computeAct_positions( dOTF, verbose)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    verbose = false;
end

[X,Y] = meshgrid(1:1:256);
R = sqrt((X-129).^2+ (Y-103).^2);
mask1 = double(R<=30);
R2 = sqrt((X-129).^2+ (Y-155).^2);
mask2 = double(R2<= 35);
overlap = mask1 .* mask2;
mask = mask1 - overlap;

load calibrated_Testbed_Acts.mat


k = (2*pi) / AOField.HeNe_Laser;
OPL = uwrap(angle(dOTF),'unwt') / k;
OPL = OPL .* mask;
OPL = OPL * 1e6;
%OPL = angle(dOTF_) / k;

for n = 1:23*23
    x = calibrated_BMC_act_locations{n}(1);
    y = calibrated_BMC_act_locations{n}(2);
   piston_area = OPL(y-1:y+1,x-1:x+1);
%     piston_area = OPL(y,x);
   Ppos(n,1) = mean(mean(piston_area(abs(piston_area)>0)));
    if isnan(Ppos(n,1))
        Ppos(n,1) = 0;
    end
end

mirror_shape = reshape(Ppos,[23,23]);
mirror_shape = padarray(mirror_shape,[5,5]);
mirror_shape = mirror_shape(1:end-1,1:end-1);
mirror_shape = circshift(mirror_shape,[-1,0]);
mirror_shape = rot90(mirror_shape,-1);
if verbose == true
    figure;
    subplot(1,2,1);
    imagesc(mirror_shape); axis xy; sqar;
    title('Computed Mirror Shape');
    
    subplot(1,2,2);
    imagesc(OPL);axis xy; sqar;
    title('OPL from dOTF');
end
end


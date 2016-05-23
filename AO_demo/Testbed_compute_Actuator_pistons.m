load('Data_Cubes_cal.mat')
dOTF = PSF_CUBE_cal.OTF - PSF_poked_CUBE_cal.OTF;
Testpattern = fitsread('DISK.fits');
Testpattern(17,7) = Testpattern(17,7) + (AOField.HeNe_Laser*1e6)/4;
[X,Y] = meshgrid(1:1:256);
R = sqrt((X-129).^2+ (Y-103).^2);
mask1 = double(R<=30);
R2 = sqrt((X-129).^2+ (Y-155).^2);
mask2 = double(R2<= 35);
overlap = mask1 .* mask2;
mask = mask1 - overlap;

load calibrated_Testbed_Acts.mat
% leftpix = [106,103];
% rightpix = [151,103];
% centerpix = [129,103];
% pokepix = [129,129];
% imagesc(rot90(Testpattern));axis xy; sqar
% xl = (linspace(leftpix(1),centerpix(1),10));
% spacing = xl(2)-xl(1);
% xl = round(linspace(leftpix(1) - 2*spacing,centerpix(1),12));
% xr = (linspace(centerpix(1),rightpix(1),10));
% spacing = xr(2)-xr(1);
% xr = round(linspace(centerpix(1),rightpix(1) + 2*spacing,12));
% y = leftpix(2) .* ones(1,12);
% %figure;plotComplex(dOTF,6)
% %axis xy
% %hold on; plot(xl,y,'r*');hold off
% %hold on; plot(xr,y,'r*');hold off
% yu = (linspace(centerpix(2),pokepix(2),11));
% spacing = yu(2)-yu(1);
% yu = round(linspace(centerpix(2),pokepix(2)+spacing,12));
% xu = centerpix(1) * ones(12,1);
% yd = (linspace(80,centerpix(2),10));
% spacing = yd(2)-yd(1);
% yd = round(linspace(80-2*spacing,centerpix(2),12));
% xd = centerpix(1) * ones(12,1);
% %hold on; plot(xd,yd,'r*');hold off
% %hold on; plot(xu,yu,'r*');hold off
% 
% 
% xr = xr(2:end);
% x = horzcat(xl,xr);
% yu = yu(2:end);
% y = horzcat(yd,yu);
% 
% [X,Y]=meshgrid(x,y);
% %R = sqrt((X-centerpix(1)).^2 + (Y-centerpix(2)).^2);
% %disk = double(R<=36);
% %X = X.*disk;
% %Y = Y.*disk;
% %plotComplex(dOTF,6)
% %axis xy
% %hold on; plot(X,Y,'g*'); hold off;
% 
% calibrated_BMC_act_locations = cell(23*23,1);
% 
% 
% counter = 1;
% for n = 1:23
% for m = 1:23
% calibrated_BMC_act_locations{counter} = [X(m,n),Y(m,n)];
% counter = counter + 1;
% end
% end

k = (2*pi) / AOField.HeNe_Laser;
OPL = uwrap(angle(dOTF),'unwt') / k;
OPL = OPL .* mask;
OPL = OPL * 1e6;
%OPL = angle(dOTF_) / k;

for n = 1:23*23
    x = calibrated_BMC_act_locations{n}(1);
    y = calibrated_BMC_act_locations{n}(2);
 %   piston_area = OPL(y-2:y+2,x-2:x+2);
    piston_area = OPL(y,x);
   Ppos(n,1) = mean(mean(piston_area(abs(piston_area)>0)));
    if isnan(Ppos(n,1))
        Ppos(n,1) = 0;
    end
end

mirror_shape = reshape(Ppos,[23,23]);
mirror_shape = padarray(mirror_shape,[5,5]);
mirror_shape = mirror_shape(1:end-1,1:end-1);
mirror_shape = circshift(mirror_shape,[-1,0]);

figure;
subplot(1,3,1);
imagesc(mirror_shape); axis xy; sqar;
title('Computed Mirror Shape');

subplot(1,3,2);
imagesc(rot90(Testpattern));axis xy; sqar;
title('Input Mirror Shape');

subplot(1,3,3);
imagesc(OPL);axis xy; sqar;
title('OPL from dOTF');
function [ dOTF, CUBE1, CUBE2] = ExtractdOTF( CUBE1, CUBE2 )
%[ dOTF ] = ExtractdOTF( CUBE1, CUBE2 )
%   Compute the dOTF from the combined Dark Corrected PSFs in the data
%   cubes created using Testbed_run_cam and ExtractAverageCUBEPSFS. This
%   will append the properties PSF_centered_and_cropped and OTF to the data
%   cubes.


PSF = CUBE1.PSF_dark_corr;
PSF_poked = CUBE2.PSF_dark_corr;

% Getting Central Pixel of PSFs
figure;
imagesc(PSF_poked);
sqar;
centerpoint1 = pickPoint(1);

% Cenering PSFs
PSF_poked = circshift(PSF_poked,1-centerpoint1);
PSF_poked = fftshift(PSF_poked);
PSF = circshift(PSF,1-centerpoint1);
PSF = fftshift(PSF);

% Crop the PSFs to 256x256 Pixels about Center Point
imagesc(PSF_poked);
sqar;
centerpoint = pickPoint(1);


PSF_poked = PSF_poked(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
PSF_poked = PSF_poked(1:end-1,1:end-1);
PSF = PSF(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
PSF = PSF(1:end-1,1:end-1);

CUBE1.PSF_centered_and_cropped = PSF;
CUBE2.PSF_centered_and_cropped = PSF_poked;

% Shift PSFs to Corners
centerpoint2 = [129,129];
PSF_poked = circshift(PSF_poked,1-centerpoint2);
PSF = circshift(PSF,1-centerpoint2);

% Compute the OTFs
OTF = fftshift(fft2(PSF_poked));
OTF(centerpoint2(1),centerpoint2(2)) = 0;
OTF_poked = fftshift(fft2(PSF));
OTF_poked(centerpoint2(1),centerpoint2(2)) = 0;

CUBE1.OTF = OTF;
CUBE2.OTF = OTF_poked;


% Compute dOTF and Store to Output
dOTF = OTF - OTF_poked;



end


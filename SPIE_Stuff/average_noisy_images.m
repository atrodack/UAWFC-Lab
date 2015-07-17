function [ Noisy_PSF1, Noisy_PSF2 ] = average_noisy_images( PSF1, PSF2, N0, Noise_Parameters )
%[ Noisy_PSF1, Noisy_PSF2 ] = average_noisy_images( PSF1, PSF2, Noise_Parameters )
% N0 = the total number of photons in the image
% Noise_Parameter Cells
% num_images = the number of images to average together
% ShotNoise = t/f flag to use Shot Noise
% ReadNoise = scale of read noise
% DarkCurrent = darkCurrent scale of detector * integration time
% N0 = number of photons in entire image

if ~iscell(Noise_Parameters)
    error('Noise_Parameters needs to be a cell array');
else
    num_images = Noise_Parameters{1};
    ShotNoise = Noise_Parameters{2};
    ReadNoise = Noise_Parameters{3};
    DarkCurrent = Noise_Parameters{4};
end

PSF1_Sum = 0;
PSF2_Sum = 0;

for n = 1:num_images
    Noise_PSF1 = addNoise(PSF1,N0,ShotNoise,ReadNoise,DarkCurrent);
    Noise_PSF2 = addNoise(PSF2,N0,ShotNoise,ReadNoise,DarkCurrent);
    PSF1_Sum = PSF1_Sum + Noise_PSF1;
    PSF2_Sum = PSF2_Sum + Noise_PSF2;
end

Noisy_PSF1 = PSF1_Sum / num_images;
Noisy_PSF2 = PSF2_Sum / num_images;


end


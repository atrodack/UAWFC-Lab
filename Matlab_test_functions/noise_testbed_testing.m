clear all;
clc;
close all;

%% Load in Images
load PSFs_from_Sim.mat;

%% Add Noise of the Same statistics to n number of pictures and sum them
<<<<<<< HEAD
num_pics = 10;
=======
num_pics = 30;
>>>>>>> e53c68d940507af119db50ab3f372eff9264ab82
PSF0_Sum = 0;
PSF1_Sum = 0;

for n = 1:num_pics
    Noisy_PSF0 = addNoise(PSF0,Field,true,1,1);
    Noisy_PSF1 = addNoise(PSF1,Field,true,1,1);
    PSF0_Sum = PSF0_Sum + Noisy_PSF0;
    PSF1_Sum = PSF1_Sum + Noisy_PSF1;
    fprintf('Picture Number %d Complete\n',n);
end

PSF0_Sum = PSF0_Sum / num_pics;
PSF1_Sum = PSF1_Sum / num_pics;

%% Do the dOTF
OTF0_Noise = fftshift(fft2(fftshift(PSF0_Sum/num_pics)));
OTF1_Noise = fftshift(fft2(fftshift(PSF1_Sum/num_pics)));


dOTF_Noise = OTF0_Noise - OTF1_Noise;
mag_Noise = abs(dOTF_Noise);
phase_Noise = angle(dOTF_Noise);


figure(1)
subplot(1,2,1)
imagesc(mag_Noise);
sqar;
colormap(gray);
title('Noisy dOTF Amplitude');

subplot(1,2,2)
imagesc(phase_Noise);
sqar;
colormap(gray);
title('Noisy dOTF Phase');



<<<<<<< HEAD
%% Generate a Mask Using Edge Detection and Image Processing

=======
% %% Generate a Mask Using Edge Detection and Image Processing
% 
>>>>>>> e53c68d940507af119db50ab3f372eff9264ab82
% [~,threshold] = edge(mag_Noise,'prewitt');
% BWs = edge(mag_Noise,'prewitt',threshold);
% se90 = strel('line',4,90);
% se0 = strel('line',4,0);
% BWsdil = imdilate(BWs,[se90 se0]);
% BWfill = imfill(BWsdil, 'holes');
% BWnobord = imclearborder(BWfill,4);
% seD = strel('disk',1,0);
% BWfinal = imerode(BWnobord,seD);
% BWfinal = imerode(BWfinal,seD);
% 
% 
% 
% 
% figure(2);
% subplot(1,2,1)
% imagesc(BWfinal.*mag_Noise);
% sqar;
% colormap(gray);
% title(sprintf('Noisy dOTF Amplitude Masked\nby\nAutogenerated Mask'));
% 
% subplot(1,2,2);
% imagesc(BWfinal.*phase_Noise);
% sqar;
% colormap(gray);
% title(sprintf('Noisy dOTF Phase Masked\nby\nAutogenerated Mask'));
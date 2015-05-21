clear all;
close all;
clc;



load Testbed_PSFs;

OTF_Finger = fftshift(fft2(img_Finger));
% OTF_Finger(241,319:323) = 0;
OTF_Finger(241,321) = 0;
OTF_No_Finger = fftshift(fft2(img_No_Finger));
% OTF_No_Finger(241,319:323) = 0;
OTF_No_Finger(241,321) = 0;

dOTF = OTF_Finger - OTF_No_Finger;
mag = abs(dOTF);
phase = angle(dOTF);

for binning = 1:10
    binning
    dOTF_binned = downsampleCCD(dOTF,binning,binning)./(binning^2);
    mag = abs(dOTF_binned);
    phase = angle(dOTF_binned);
    
    figure(1)
    subplot(1,3,1);
    plotComplex(dOTF_binned,2); drawnow;
    sqar;
    axis xy;
    
    subplot(1,3,2);
    imagesc(mag);
    sqar;
    axis xy;
    
    subplot(1,3,3);
    imagesc(phase);
    sqar;
    axis xy;
    colormap(gray);
    input 'more?';
end



% dOTF_binned = downsampleCCD(dOTF/16,4,4);
% mag = abs(dOTF_binned);
% phase = angle(dOTF_binned);
% 
% uwrapphase = uwrap(phase,'flyn');
% 
% 
% figure(1);
% 
% 
% figure;
% imagesc(uwrapphase);
% axis xy;
% sqar;
% colormap(gray);

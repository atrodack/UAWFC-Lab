clear all;
clc;
close all;


load PSFs_from_Sim.mat;

OTF0 = fftshift(fft2(fftshift(PSF0)));
OTF1 = fftshift(fft2(fftshift(PSF1)));

dOTF = OTF0 - OTF1;
mag = abs(dOTF);
phase = angle(dOTF);

[~,threshold] = edge(mag,'sobel');
BWs = edge(mag,'sobel',threshold);

Gx = [-1,1];
Gy = Gx';
Ix = conv2(mag,Gx,'same');
Iy = conv2(mag,Gy,'same');
Grad_mag = sqrt(Ix.^2 + Iy.^2);
Grad_mag = Grad_mag/max(Grad_mag(:));

se90 = strel('line',4,90);
se0 = strel('line',4,0);
BWsdil = imdilate(Grad_mag,[se90 se0]);
BWfill = imfill(BWsdil, 'holes');
BWnobord = imclearborder(BWfill,4);
seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);




figure;
imagesc(BWfinal);
sqar;


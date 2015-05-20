clear all;
close all;
clc;



Num_Folders = 2;
Num_Files_per_Folder = 100;

% varargin{1} = '/home/alex/Desktop/Git_Repos/UAWFC/Real_Data/201558_Batch1_bp10_PSFWithFinger';
% varargin{2} = '/home/alex/Desktop/Git_Repos/UAWFC/Real_Data/201558_Batch2_bp10_PSFWithoutFinger';
varargin{1} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch1_bp10_PSFWithFinger/';
varargin{2} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch2_bp10_PSFWithoutFinger/';
varargin{3} = 'RAW_scienceIM_frame_';
varargin{4} = 'RAW_scienceIM_frame_';

imgs = BatchRead(Num_Folders,Num_Files_per_Folder,false,varargin);
close;

img_Finger = AddImages(imgs{1});
img_No_Finger = AddImages(imgs{2});

imagesc(img_Finger);
sqar;


center_PSF = pickPoint(1);
img_Finger = circshift(img_Finger,1-center_PSF);
img_No_Finger = circshift(img_No_Finger,1-center_PSF);

OTF_Finger = fftshift(fft2(img_Finger));
OTF_Finger(241,321) = 0;
OTF_No_Finger = fftshift(fft2(img_No_Finger));
OTF_No_Finger(241,321) = 0;

dOTF = OTF_Finger - OTF_No_Finger;
mag = abs(dOTF);
phase = angle(dOTF);


uwrapphase = uwrap(phase,'flyn');


close all;
figure(1);
subplot(1,2,1);
imagesc(mag);
sqar;
axis xy;

subplot(1,2,2);
imagesc(uwrapphase);
sqar;
axis xy;
colormap(jet);
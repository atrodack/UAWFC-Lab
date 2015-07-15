% Read in data from Real_Data folder
clear all;
clc;
close all;


Num_Folders = 2;
Num_Files_per_Folder = 10;

% varargin{1} = '/home/alex/Desktop/Git_Repos/UAWFC/Real_Data/201558_Batch1_bp10_PSFWithFinger';
% varargin{2} = '/home/alex/Desktop/Git_Repos/UAWFC/Real_Data/201558_Batch2_bp10_PSFWithoutFinger';
varargin{1} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch1_bp10_PSFWithFinger/';
varargin{2} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch2_bp10_PSFWithoutFinger/';
varargin{3} = 'RAW_scienceIM_frame_';
varargin{4} = 'RAW_scienceIM_frame_';

imgs = BatchRead(Num_Folders,Num_Files_per_Folder,false,varargin);
close;

img_Finger = AddImages(imgs{1});
% img_Finger(img_Finger<315) = 0;
img_No_Finger = AddImages(imgs{2});
% img_No_Finger(img_No_Finger<315) = 0;


figure(1);
imagesc(img_Finger);
sqar;
title('Summed Finger Image');

figure(2);
imagesc(img_No_Finger);
sqar;
title('Summed No Finger Image');


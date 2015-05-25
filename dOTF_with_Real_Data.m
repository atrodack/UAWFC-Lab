clear all;
close all;
clc;


Num_Folders = 2;
Num_files_per_folder = 100;
varargin{1} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch1_bp10_PSFWithFinger/';
varargin{3} = 'RAW_scienceIM_frame_';
varargin{2} = '/home/alex/Desktop/UAWFC/Real_Data/2015520_Batch2_bp10_PSFWithoutFinger/';
varargin{4} = 'RAW_scienceIM_frame_';

images = BatchRead(Num_Folders,Num_files_per_folder, false, varargin);

img_Finger = images{1};
img_Finger = AddImages(img_Finger);
img_Finger = img_Finger(34:290,164:420);
img_Finger = img_Finger(1:end-1,1:end-1);
% img_Finger = img_Finger(114:210,228:356);
% img_Finger = padarray(img_Finger,[(210-114)/2,(356-228)/2],'both');
img_No_Finger = images{2};
img_No_Finger = AddImages(img_No_Finger);
img_No_Finger = img_No_Finger(34:290,164:420);
img_No_Finger = img_No_Finger(1:end-1,1:end-1);
% img_No_Finger = img_No_Finger(114:210,228:356);
% img_No_Finger = padarray(img_No_Finger,[(210-114)/2,(356-228)/2],'both');

imagesc(img_Finger);
sqar;
pt = pickPoint(1);

img_Finger = circshift(img_Finger,1-pt);
img_No_Finger = circshift(img_No_Finger,1-pt);


OTF_Finger = fftshift(fft2(img_Finger));
% OTF_Finger(241,319:323) = 0;
% OTF_Finger(241,321) = 0;
OTF_Finger(pt(1),pt(2)) = 0;
OTF_No_Finger = fftshift(fft2(img_No_Finger));
% OTF_No_Finger(241,319:323) = 0;
% OTF_No_Finger(241,321) = 0;
OTF_No_Finger(pt(1),pt(2)) = 0;

dOTF = OTF_Finger - OTF_No_Finger;
mag = abs(dOTF);
phase = angle(dOTF);


clf;
imagesc(phase);
pt2 = pickPoint(1);
[X,Y] = mkImageCoords(dOTF,1,pt2);
R = sqrt(X.^2+Y.^2);

%%
for alpha = linspace(-1,0,1000)
    defocus = exp(alpha*1i.*R);
    newphase = angle(defocus.*dOTF);
%     unwrappedphase = uwrap(newphase,'unwt');
%     imagesc(unwrappedphase);
imagesc(newphase);
    title(sprintf('alpha = %0.3f ',alpha));
    drawnow;
end



% for binning = 1:10
%     binning
%     dOTF_binned = downsampleCCD(dOTF,binning,binning)./(binning^2);
%     mag = abs(dOTF_binned);
%     phase = angle(dOTF_binned);
%     
%     figure(1)
%     subplot(1,3,1);
%     plotComplex(dOTF_binned,2); drawnow;
%     sqar;
%     axis xy;
%     
%     subplot(1,3,2);
%     imagesc(mag);
%     sqar;
%     axis xy;
%     
%     subplot(1,3,3);
%     imagesc(phase);
%     sqar;
%     axis xy;
%     colormap(gray);
%     input 'more?';
% end


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

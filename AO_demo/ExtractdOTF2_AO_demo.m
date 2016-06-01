function [ dOTF, CUBE1, CUBE2] = ExtractdOTF2_AO_demo( CUBE1, CUBE2, version)
%[ dOTF ] = ExtractdOTF( CUBE1, CUBE2 )
%   Compute the dOTF from the combined Dark Corrected PSFs in the data
%   cubes created using Testbed_run_cam and ExtractAverageCUBEPSFS. This
%   will append the properties PSF_centered_and_cropped and OTF to the data
%   cubes.

% Make old scripts work
if nargin < 3
    version = 1;
end

if version == 1
    %% OLD WAY
    PSF = CUBE1.PSF_dark_corr;
    PSF_poked = CUBE2.PSF_dark_corr;
    
%   HeNe (needs recalibration after movement of flip mirror/OD filter)
%   centerpoint = [272,305]; %[y,x] in ij
    
%   SuperK (needs recalibration after input of OD filter)
    centerpoint = [254,316]; %no OD filter
    
    % Cenering PSFs and cropping to 256x256 (remove tilt)
    PSF_poked = PSF_poked(centerpoint(1)-127:centerpoint(1)+128,centerpoint(2)-127:centerpoint(2)+128);
    PSF = PSF(centerpoint(1)-127:centerpoint(1)+128,centerpoint(2)-127:centerpoint(2)+128);
    
    CUBE1.PSF_centered_and_cropped = PSF;
    CUBE2.PSF_centered_and_cropped = PSF_poked;
    
    % Shift PSFs to Corners (fftshift by hand to make sure pixel allignment is correct)
    centerpoint2 = [129,129]; %for 256x256
    PSF_poked = circshift(PSF_poked,1-centerpoint2);
    PSF = circshift(PSF,1-centerpoint2);
    
    % Compute the OTFs and zero central pixel
    OTF = fftshift(fft2(PSF_poked));
    OTF(centerpoint2(1),centerpoint2(2)) = 0;
    OTF_poked = fftshift(fft2(PSF));
    OTF_poked(centerpoint2(1),centerpoint2(2)) = 0;
    
    CUBE1.OTF = OTF;
    CUBE2.OTF = OTF_poked;
    
    
    % Compute dOTF and Store to Output
    dOTF = OTF - OTF_poked;
    
elseif version == 2
    %% NEW WAY
    [CUBE1, CUBE2, ~] = ExamineCUBEs(CUBE1, CUBE2, 2); % Compute OTF Cubes
    [CUBE1, CUBE2, ~] = ExamineCUBEs(CUBE1, CUBE2, 4); % Compute average OTF and PSF
    [CUBE1, CUBE2, dOTF] = ExamineCUBEs(CUBE1, CUBE2, 5); % Compute dOTF
%     [CUBE1, CUBE2, ~] = ExamineCUBEs(CUBE1, CUBE2, 6); % align the OTF stacks
%     [CUBE1, CUBE2, ~] = ExamineCUBEs(CUBE1, CUBE2, 8); % Compute average of shifted OTF stacks
%     [CUBE1, CUBE2, dOTF] = ExamineCUBEs(CUBE1, CUBE2, 11); % Scan Power to equalized and compute dOTF
    
    
    %% NOTES for DATA CUBE
    note1 = sprintf('PSFs is the raw data for PSFs from the camera');
    note2 = sprintf('DARKS is the raw data for Dark images from the camera. If none are taken, they are just matrices of zeroes');
    note3 = sprintf('Camera_Data is a struct holding info returned while camera is taking images');
    note4 = sprintf('PSF_dark_corr is the dark corrected, average PSF (size 480x640)');
    note5 = sprintf('cropped_PSFs is the data for PSFs centered at the calibrated position and cropped to size 256x256');
    note6 = sprintf('OTFs is the data for the OTFs computed for cropped_PSFs');
    note7 = sprintf('PSF is the average of cropped_PSFs frames');
    note8 = sprintf('OTF is the average of OTFs frames');
%     note9 = sprintf('TT is the Tip/Tilt outputs of vernierAlignOTF');
%     note10a = sprintf('OTF_ is the average of the vernier aligned OTFs (multiplied by Scale_Factor)');
%     note10b = sprintf('OTF_ is the average of the vernier aligned OTFs');
%     note11 = sprintf('Scale_Factor is the computed power equalizing factor for CUBE1 to CUBE2');
    
    NOTE1 = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',note1,note2,note3,note4,note5,note6,note7,note8);
    NOTE2 = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',note1,note2,note3,note4,note5,note6,note7,note8);
    
    CUBE1.NOTE = NOTE1;
    CUBE2.NOTE = NOTE2;
end %if version


end %function


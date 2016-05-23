function [CUBE] = ExtractAverageCUBEPSFs_AO_demo( CUBE )
%[CUBE] = ExtractAverageCUBEPSFs( CUBE )
%   Takes a data cube from Testbed_run_cam and computes the combined PSF
%   and combined Dark Corrected PSF, and appends them to the data structure

nframes = size(CUBE.PSFs,4);


%% Uncorrected Average Image
IMG = 0;

for n = 1:nframes
    IMG = double(IMG) + double(CUBE.PSFs(:,:,1,n));
end

IMG = IMG / nframes;

%% Dark Corrected Average Image
IMG_CORR = 0;

for n = 1:nframes
    IMG_CORR = double(IMG_CORR) + double(CUBE.PSFs(:,:,1,n)) - double(CUBE.DARKS(:,:,1,n));
end

IMG_CORR = IMG_CORR / nframes;



CUBE.PSF = IMG;
CUBE.PSF_dark_corr = IMG_CORR;

end


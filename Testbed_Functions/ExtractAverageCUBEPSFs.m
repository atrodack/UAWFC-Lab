function [CUBE] = ExtractAverageCUBEPSFs( CUBE )
%[CUBE] = ExtractAverageCUBEPSFs( CUBE )
%   Takes a data cube from Testbed_run_cam and computes the combined PSF
%   and combined Dark Corrected PSF, and appends them to the data structure

nframes = size(CUBE.PSFs,3);


%% Uncorrected Average Image
IMG = 0;

for n = 1:nframes
    IMG = IMG + CUBE.PSFs(:,:,n);
end

IMG = IMG / nframes;

%% Dark Corrected Average Image
IMG_CORR = 0;

for n = 1:nframes
    IMG_CORR = IMG_CORR + CUBE.PSFs(:,:,n) - CUBE.DARKS(:,:,n) - double(CUBE.BACKGROUND);
end

IMG_CORR = IMG_CORR / nframes;



CUBE.PSF = IMG;
CUBE.PSF_dark_corr = IMG_CORR;

end


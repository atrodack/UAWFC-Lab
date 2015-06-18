function Noisy_PSF = addNoise(PSF, Field, ShotNoise, ReadNoise, DarkCurrent )
%Noisy_PSF = addNoise(PSF, Field, ShotNoise, ReadNoise, DarkCurrent )
%
% Inputs
% PSF = a Point Spread Function
% Field = the grid of an AOField Object used to get PSF
% ShotNoise = t/f flag to use Shot Noise
% ReadNoise = scale of read noise
% DarkCurrent = darkCurrent scale of detector * integration time

if nargin < 4
    if ShotNoise == true
        Sum0 = sum(PSF(:));
        N0 = sum(sum(abs(Field)));
        Noisy_PSF = PSF*(1e-12*N0/Sum0);
        Noisy_PSF = 1e12*imnoise(Noisy_PSF,'poisson');
    end
    
elseif nargin < 5
    if ShotNoise == true
        Sum0 = sum(PSF(:));
        N0 = sum(sum(abs(Field)));
        Noisy_PSF = PSF*(1e-12*N0/Sum0);
        Noisy_PSF = 1e12*imnoise(Noisy_PSF,'poisson');
    end
    Noisy_PSF = Noisy_PSF + randn(size(PSF)) * ReadNoise;
    
elseif nargin < 6
    if ShotNoise == true
        Sum0 = sum(PSF(:));
        N0 = sum(sum(abs(Field)));
        Noisy_PSF = PSF*(1e-12*N0/Sum0);
        Noisy_PSF = 1e12*imnoise(Noisy_PSF,'poisson');
    end
    Noisy_PSF = Noisy_PSF + randn(size(PSF)) * ReadNoise;
    Noisy_PSF = Noisy_PSF + randn(size(PSF)) * DarkCurrent;
end

    











end
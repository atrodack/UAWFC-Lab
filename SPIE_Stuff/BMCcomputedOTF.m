function [ dotf, psf1, psf2, otf1, otf2, dotfd] = BMCcomputedOTF( DM, pokeact, actuators_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2, DECONVOLVE )
%[ dOTF, PSF1, PSF2, OTF1, OTF2 ] = BMCcomputedOTF( DM, pokeact, actuators_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2, DECONVOLVE )
%   
load dOTF_act_698_mask.mat

if nargin < 7
    Phasescreen1 = 1;
    Phasescreen2 = 1;
    DECONVOLVE = false;
elseif nargin < 8
    Phasescreen2 = 1;
    DECONVOLVE = false;
elseif nargin < 9
    DECONVOLVE = false;
end

if ~iscell(Noise_Parameters)
    if islogical(Noise_Parameters)
        usenoise = Noise_Parameters;
        Noise_Parameters = cell(5,1);
        Noise_Parameters{5} = usenoise;
    else
        error('Noise_Parameters must be a cell array');
    end
else
    if isempty(Noise_Parameters{5})
        Noise_Parameters{5} = false; %if flag is unset, don't use noise
        N0 = Noise_Parameters{6};
    else
        N0 = Noise_Parameters{6};
    end
end

lambda = Field.lambda;



F1 = Field.copy;
F1.name = 'BMC Field 1';

F2 = Field.copy;
F2.name = 'BMC Field 2';

% display('Initializing DM dOTF Positions');

Ppos_flat = actuators_in;
Ppos_poked = actuators_in;
Ppos_poked(pokeact) = Ppos_poked(pokeact) + (lambda/4);


% display('Flattening DM');

DM.setActs(Ppos_flat);
DM.touch;
DM.render;


F1.planewave * Phasescreen1 * Phasescreen2 * A * DM;
grid1 = F1.grid;
psf1 = abs(fftshift(fft2(fftshift(grid1)))).^2;

% PSF1 = F1.mkPSF(FOV,PLATE_SCALE);

% fprintf('Poking Segment %d \n',pokeseg);

DM.setActs(Ppos_poked);
DM.touch;
DM.render;

% display('Making PSF 2');

F2.planewave * Phasescreen1 * Phasescreen2 * A * DM;
grid2 = F2.grid;
psf2 = abs(fftshift(fft2(fftshift(grid2)))).^2;

% PSF2 = F2.mkPSF(FOV,PLATE_SCALE);

% if Noise_Parameters{5} == true
%     [ PSF1, PSF2 ] = average_noisy_images( PSF1, PSF2, Field.grid, Noise_Parameters );
% end
    
if Noise_Parameters{5} == true
    [ psf1, psf2 ] = average_noisy_images( psf1, psf2, N0, Noise_Parameters );
end

% F1.touch; F2.touch;
% F1.grid(PSF1); F2.grid(PSF2);

% display('Making OTFs');
% OTF1 = F1.mkOTF2(FoV,new_spacing(1));
% OTF2 = F2.mkOTF2(FoV,new_spacing(1));

otf1 = fftshift(fft2(fftshift(psf1))) .* (1 / F1.dx)^2;
otf1(1025,1025) = 0;
otf2 = fftshift(fft2(fftshift(psf2))) .* (1 / F2.dx)^2;
otf2(1025,1025) = 0;



% F1.touch; F2.touch;

% display('Computing dOTF');

% dOTF = OTF1 - OTF2;
dotf = otf1 - otf2;



%% Deconvolution
if DECONVOLVE == true
    fdiff = grid1 - grid2;
    % figure;
    % plotComplex(fdiff,2);
    load ALEX_P.mat;
    FDIFF = fftshift(fft2(circshift(conj(fdiff),1-ALEX_P)));
    DOTF = fftshift(fft2(fftshift(dotf)));

    gamma = 1.5e5;
    Deconv = Wiener(DOTF,conj(FDIFF),10,gamma);
    
    dotfd = Deconv;
else
    dotfd = [];

end



DM.setActs(Ppos_flat);
DM.touch;
DM.render;
end



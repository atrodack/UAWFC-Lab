function [ dOTF, PSF1, PSF2, OTF1, OTF2, dotfd ] = IrisAOcomputedOTF_new( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2, NECO, DECONVOLVE )
%[ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2 )
%   

if nargin < 7
    Phasescreen1 = 1;
    Phasescreen2 = 1;
    NECO = false;
elseif nargin < 8
    Phasescreen2 = 1;
    NECO = false;
elseif nargin < 9
    NECO = false;
end

if ~iscell(Noise_Parameters)
    error('Noise_Parameters must be a cell array');
else
    if isempty(Noise_Parameters{5})
        Noise_Parameters{5} = false; %if flag is unset, don't use noise
        N0 = Noise_Parameters{6};
    else
        N0 = Noise_Parameters{6};
    end
end


new_spacing = DM.spacing;
dx = new_spacing(1);

F1 = Field.copy;
F1.name = 'IrisAO Field 1';

F2 = Field.copy;
F2.name = 'IrisAO Field 2';

% display('Initializing DM dOTF Positions');

PTTpos_flat = PTTpos_in;
PTT_flat = mapSegments(PTTpos_flat);
PTTpos_poked = PTTpos_in;



%
if NECO == false
    PTTpos_poked(pokeseg,1) = (F1.lambda) / 4;
else
    PTTpos_poked(pokeseg,3) = 0.1e-3;
    PTTpos_poked(pokeseg,1) = (F1.lambda) / 4;
end
%



PTT_poked = mapSegments(PTTpos_poked);

% display('Flattening DM');
DM.setIrisAO(PTT_flat);

F1.planewave * Phasescreen1 * Phasescreen2 * A * DM;
grid1 = F1.grid;
PSF1 = abs(fftshift(fft2(fftshift(grid1)))*(dx^2)).^2;

DM.setIrisAO(PTT_poked);

F2.planewave * Phasescreen1 * Phasescreen2 * A * DM;
grid2 = F2.grid;
PSF2 = abs(fftshift(fft2(fftshift(grid2)))*(dx^2)).^2;

if Noise_Parameters{5} == true
    [ PSF1,PSF2 ] = average_noisy_images( PSF1, PSF2, N0, Noise_Parameters );
end


OTF1 = fftshift(fft2(fftshift(PSF1)))*((1/dx)^2);
OTF2 = fftshift(fft2(fftshift(PSF2))*((1/dx)^2));

dOTF = OTF1 - OTF2;


%% Deconvolution
if DECONVOLVE == true
    fdiff = grid1 - grid2;
%     figure;
%     plotComplex(fdiff,2);
%     ALEX_P = pickPoint;
%     save('ALEX_P_IRISAO','ALEX_P');
    load ALEX_P_IRISAO.mat;
    FDIFF = fftshift(fft2(circshift(conj(fdiff),1-ALEX_P)));
    DOTF = fftshift(fft2(fftshift(dOTF)));

    gamma = 8.5e3;
    Deconv = Wiener(DOTF,conj(FDIFF),10,gamma);
    
    dotfd = Deconv;
else
    dotfd = [];

end

end
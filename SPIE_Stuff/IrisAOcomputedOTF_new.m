function [ dOTF, PSF1, PSF2, OTF1, OTF2, dotfd ] = IrisAOcomputedOTF_new( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2, NECO, DECONVOLVE,bandpass )
%[ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2 )
%   
dOTF = [];
OTF1 = [];
OTF2 = [];

if nargin < 7
    Phasescreen1 = 1;
    Phasescreen2 = 1;
    NECO = false;
    DECONVOLVE = false;
    banpass = Field.lambda;
elseif nargin < 8
    Phasescreen2 = 1;
    NECO = false;
    DECONVOLVE = false;
    bandpass = Field.lambda;
elseif nargin < 9
    NECO = false;
    DECONVOLVE = false;
    bandpass = Field.lambda;
elseif nargin < 10
    bandpass = Field.lambda;
    DECONVOLVE = false;
elseif nargin < 11
    bandpass = Field.lambda;
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
% PSFBAND1 = 0;
% PSFBAND2 = 0;
PSFBAND1 = cell(1,0);
PSFBAND2 = cell(1,0);
dOTFBAND = 0;
FoV = F1.FOV;
n = 1;
for lambda_passband = 1:length(bandpass)
    % display('Flattening DM');
    DM.setIrisAO(PTT_flat);
    F1.lambda = bandpass(lambda_passband);
    fprintf('.');
    F1.planewave * Phasescreen1 * Phasescreen2 * A * DM;
    grid1 = F1.grid;
    psf1 = (fftshift(fft2(fftshift(grid1)))*(dx^2));
    
    k = 2*pi/F1.lambda;
    [kx,ky] = kcoords(F1);
    
    % This is the size of the FFT pixels...
    dTH = F1.dk/F1.k*206265;
    thx_ = kx/k*206265;
    thy_ = ky/k*206265;
    
    dth = median(diff(thx_));
    
    
    % this is get a slightly larger region to interpolate into.
    % No more extrapolation NaNs!
    SELx = abs(thx_)<=(FoV+2*dTH(1));
    SELy = abs(thy_)<=(FoV+2*dTH(2));
    thx_ = thx_(SELx);
    thy_ = thy_(SELy);
    [THX_,THY_] = meshgrid(thx_,thy_);
    
    Nth = ceil(FoV/dth);
    
    thx = (-Nth:Nth)*dth;
    thy = thx;
    [THX,THY] = meshgrid(thx,thy);
    
    psf1 = psf1(SELx,SELy);
    
    HALO = qinterp2(THX_,THY_,psf1,THX,THY);
    PSF1 = abs(HALO).^2;
    PSF1(isnan(PSF1)) = 0;
    
%     PSFBAND1 = PSFBAND1 + PSF1;
    PSFBAND1{n} = PSF1;
    
    DM.setIrisAO(PTT_poked);
    F2.lambda = bandpass(lambda_passband);
    F2.planewave * Phasescreen1 * Phasescreen2 * A * DM;
    grid2 = F2.grid;
    psf2 = (fftshift(fft2(fftshift(grid2)))*(dx^2));
    
    k = 2*pi/F2.lambda;
    [kx,ky] = kcoords(F2);
    
    % This is the size of the FFT pixels...
    dTH = F2.dk/F2.k*206265;
    thx_ = kx/k*206265;
    thy_ = ky/k*206265;
    
    dth = median(diff(thx_));
    
    
    % this is get a slightly larger region to interpolate into.
    % No more extrapolation NaNs!
    SELx = abs(thx_)<=(FoV+2*dTH(1));
    SELy = abs(thy_)<=(FoV+2*dTH(2));
    thx_ = thx_(SELx);
    thy_ = thy_(SELy);
    [THX_,THY_] = meshgrid(thx_,thy_);
    
    Nth = ceil(FoV/dth);
    
    thx = (-Nth:Nth)*dth;
    thy = thx;
    [THX,THY] = meshgrid(thx,thy);
    
    psf2 = psf2(SELx,SELy);
    
    HALO = qinterp2(THX_,THY_,psf2,THX,THY);
    PSF2 = abs(HALO).^2;
    PSF2(isnan(PSF2)) = 0;
%     PSFBAND2 = PSFBAND2 + PSF2;
    PSFBAND2{n} = PSF2;

    
    if Noise_Parameters{5} == true
        [ PSF1,PSF2 ] = average_noisy_images( PSF1, PSF2, N0, Noise_Parameters );
    end
    
    
%     OTF1 = fftshift(fft2(fftshift(PSF1)))*((1/dx)^2);
%     OTF1(1025,1025) = 0;
%     OTF2 = fftshift(fft2(fftshift(PSF2))*((1/dx)^2));
%     OTF2(1025,1025) = 0;
%     dOTF = OTF1 - OTF2;
%     dOTFBAND = dOTFBAND + dOTF;
n = n+1;
end
fprintf('\n');
% PSF1 = PSFBAND1 / length(bandpass);
% PSF2 = PSFBAND2 / length(bandpass);
% dOTF = dOTFBAND / length(bandpass);
PSF1 = PSFBAND1;
PSF2 = PSFBAND2;

DM.setIrisAO(PTT_flat);

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

%     gamma = 8.5e3;
gamma = 2.5e5;
    Deconv = Wiener(DOTF,conj(FDIFF),10,gamma);
    
    dotfd = Deconv;
else
    dotfd = [];

end

end
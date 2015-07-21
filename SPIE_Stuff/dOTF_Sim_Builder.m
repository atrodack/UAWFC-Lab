% function Estimate = dOTF_Sim_Builder(Estimate)
clear all;
clc;
close all;
% Deconvolution Testing

%**************************************************************************
%                       SPIE dOTF Deconvolution Methods Simulation
%**************************************************************************

%% System Parameters
% Don't touch these if using maketestbedelements.m to construct components
global lambda k D secondary spider SPACING aa fftsize THld FOV PLATE_SCALE FoV_withIrisAO FoV_withoutIrisAO RunSIM RunTESTBED IrisAO_on BMC_on verbose_makeDM Scalloped_Field UseRealPSF coronagraph system_verbose

% Set Wavelength
% lambda = AOField.RBAND; % Red light
lambda = AOField.HeNe_Laser; %Testbed wavelength

% Compute Wavenumber
k = (2*pi) / lambda;

% Pupil Specs
D = 3.636e-3; % width of the longest row of segments for a 37 segment IrisAO mirror
% secondary = 0.3*D; % 30% of Pupil Diameter
secondary = 0;
% spider = 0.02*D; % 2% of Pupil Diameter
spider = 0;


%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing
nzerns = 4; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
fftsize = 2^10;

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs
FOV = 100*THld; % arcsecs
PLATE_SCALE = THld/4;
FoV_withIrisAO = 3.5e-3;
FoV_withoutIrisAO = 10.5e-3;

%% Run Simulation/Testbed Flags
RunSIM = true; %Run the simulation
RunTESTBED = false; %Run the testbed equipment
Plotting = true; % turn on/off plotting interesting pictures of the dOTF estimate, also dictates if debugging will print pictures
debuggery = true; % Run intermediate steps of interest for debugging!
looping = false; % change dOTF_Sim_Builder into a function which outputs pertinent data for deconvolution studies/debugging.  If looping is true, must manually uncomment the function and end lines.

%% Simulation Flags

% IrisAO Flags
IrisAO_on = true; % turns on/off IrisAO Mirror (if false, DM1 variable set to 1)
verbose_makeDM = false; % turns on/off plotting the mirror as it is constructed
Scalloped_Field = false; % turns on/off returning an AOField Object that encodes the actual surface shape of the segments
spiral = true; % turns on/off pattern used in NECO data
segment = true; % use a segment to generate the difference field.  Default to true, but if it is false, use a blocker to compute dOTF.
exactDifference = false; % choose to use the exact difference, or compute the difference from the field objects - needs to be updated to be more efficient with what to calculate and to deal with partial obscurations
Minus = true; % determines portion of the dOTF you wish to mask and deconvolve, default is the lower field (referred to as dO_Plus)
manualPick = false; % turns on/off the ability to select pixel points for recentering of the difference field and the masked portion of the dOTF
pointOverwrite = false; % only accessible when manualPick is true, in which case it will save the selected points for future use (overwrites the current point holding files)
smoothing = true; % turns on/off smoothing the mask edge applied to the dOTF

% BMC Flag
BMC_on = false; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = false; %Injects nzerns Zernike Terms
InjectRandAb = true; %if InjectAB is true, picks Zernikes "Randomly"
InjectKnownAb = false; %if InjectAB is true, picks provided Zernikes

InjectKolm = false;

% Check Aberration Flags
if InjectKolm == true
    if InjectAb == true
        y = 'y';
        Y = 'y';
        n = 'n';
        N = 'n';
        display('You have turned on both Zernikes and Kolmogorov Turbulence. Would you like to turn one off?');
        answer1 = input('y/n : ');
        if strcmpi(answer1,'y') == true
            Zern = 'Zern';
            zern = 'Zern';
            Kolm = 'Kolm';
            kolm = 'Kolm';
            display('Turn off which Aberration Type?');
            answer2 = input('Zern/Kolm : ');
            if strcmpi(answer2,'Zern') == true
                InjectAb = false;
                display('Zernikes will no longer be injected');
            elseif strcmpi(answer2,'Kolm') == true
                InjectKolm = false;
                display('Kolmogorov Turbulence will no longer be injected');
            else
                error('Incorrect Input: Must be either Zern or Kolm');
            end
            clear Zern Kolm zern kolm
        end
        clear y n Y N
    end
end

% Corrector Flag
UseDM4Correction = true;

% Noise Flags
UseNoise = false;
if UseNoise == true
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 1;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 10;
    Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = 5e6;
else
    Noise_Parameters = cell(5,1);
    %     Noise_Parameters{1} = 1;
    %     Noise_Parameters{2} = false;
    %     Noise_Parameters{3} = 0;
    %     Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
end

% Use Testbed PSF Instead of Simulated PSF
UseRealPSF = false;
if UseRealPSF == true
    InjectAb = false;
    Num_Folders = 2;
    Num_files_per_folder = 100;
    %     varargin{1} = '/home/alex/Desktop/Data/2015612_Batch1_nofilter_PSFWithoutFinger/';
    %     varargin{3} = 'RAW_scienceIM_frame_';
    %     varargin{2} = '/home/alex/Desktop/Data/2015612_Batch2_nofilter_PSFWithFinger/';
    %     varargin{4} = 'RAW_scienceIM_frame_';
    %     varargin{1} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithoutFingerDMBox/';
    %     varargin{3} = 'RAW_scienceIM_frame_';
    %     varargin{2} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithFingerDMBox/';
    %     varargin{4} = 'RAW_scienceIM_frame_';
end
% Coronagraph Flag
coronagraph = false; % turns on going through coronagraph elements

% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Make the Testbed Elements
% Or load them in...
numRings = 3;
numSeg = sum(1:numRings)*6 + 1;
load MadeTestbedElements.mat;
% maketestbedelements;
A.lambdaRef = lambda;
% Initialize exact difference field
DMDiff = makeIrisAODM(1,verbose_makeDM,Scalloped_Field,0);
DMDiff.lambdaRef = lambda;
%% Inject Aberration
if RunSIM == true
    if InjectAb == true && InjectRandAb == true
        ABER = AOScreen(A);
        n = sort((randi(4,1,nzerns)),'ascend'); %zernike mode order (from lowest to highest)
        m = zeros(1,nzerns); %zernike mode initialization
        
        %Get a correct "m" index for each "n" index
        for ii = 1:nzerns
            m_pos = -n(ii):2:n(ii);
            choice = randi(length(m_pos),1,1);
            m(ii) = m_pos(choice);
        end
        
        % Find a coefficient between -1 and 1 waves, scale it to the number of
        % zernikes added to avoid unrealistic PSFs (fairly arbitrary,
        % probably uncessary)
        if nzerns < 3
            coeffs = (2*rand(1,nzerns)-1);
        elseif nzerns <= 5
            coeffs = (((nzerns-(nzerns-2))/nzerns)) .* (2*rand(1,nzerns)-1);
        elseif nzerns > 5
            coeffs = (((nzerns-(nzerns-5))/nzerns)) .* (2*rand(1,nzerns)-1);
        end
        
        % Add the Zernikes into ABER
        ABER.zero;
        for ii = 1:nzerns
            ABER.addZernike(n(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n = n';
        m = m';
        Number_of_waves = coeffs';
        T = table(n,m,Number_of_waves);
        fprintf('\nInjected Aberrations:\n');
        disp(T);
    elseif InjectAb == true && InjectKnownAb == true
        ABER = AOScreen(A);
        %         n = [2,2,2,3,3];
        n = [1,1,4];
        %         m = [-2,0,2,-1,3];
        m = [-1,1,0];
        
        %         coeffs = 1 * randn(1,length(n));
        %         coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
        coeffs = 0.25*randn(1,length(n));
        ABER.zero;
        for ii = 1:length(n)
            ABER.addZernike(n(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n = n';
        m = m';
        Number_of_waves = coeffs';
        T = table(n,m,Number_of_waves);
        fprintf('\nInjected Aberrations:\n');
        disp(T);
        ABER.show;
        aberration = ABER.grid;
    elseif InjectAb == false
        ABER = 1;
    end
    
     if InjectKolm == true
        TURB = AOScreen(A,4e-3,lambda);
        TURB.spacing(SPACING);
        TURB.name = 'Simulated Turbulence';
        TURB.make;
        %         grid = TURB.grid;
        %         TURB.grid(grid * 30);
        wind_dir = randn(2,1);
        wind_dir = wind_dir./abs(wind_dir);
        wind_strength = randi(5,2,1);
        Wind = wind_dir .* wind_strength;
        
        r0 = TURB.r0;
        D_fit = D;
        %  expected_strehl = exp(-(1.03 * (D_fit/r0) ^ (5/3)))
        expected_strehl = exp(-(0.13 * (D_fit/r0) ^ (5/3)));
    else
        TURB = 1;
     end
end

%% IrisAO Simulation
F = AOField(fftsize);
F.FFTSize = (fftsize);
F.spacing(SPACING);
F.lambda = lambda;
F.name = 'IrisAO Field 1';

F2 = F.copy;
F2.name = 'IrisAO Field 2';

FD = F.copy;
FD.name = 'IrisAO Difference Field';

mask = F.copy;
fd.name = 'dOTF Noise Mask';

% Initialize Segment Locations
% PTTpos_flat = zeros(numSeg,3);
% PTT_flat = mapSegments(PTTpos_flat,numRings);

PTTpos_poked1 = zeros(numSeg,3);

% Generate piston spiral pattern on 1st 7 segments on DM surface
if spiral == 1
    piston_spiral = 25e-9; % 25 nm piston
    for seg = 1:7
        PTTpos_poked1(seg,1) = piston_spiral;
        piston_spiral = piston_spiral + 25e-9;
    end
end

PTTpos_poked2 = PTTpos_poked1; % copy the spiral pattern

if segment == true
    % Create the location (segment) of the difference field
    if numRings > 0
        DiffField = 3;
        for segRing = 1:numRings-1
            DiffField = DiffField + (segRing)*6 + 1;
        end
    else
        DiffField = 1;
    end
    
    % Set Difference Field Position(s)
    PTT1 = logical([0 0 1]); % turns on/off different piston, tip/tilt difference field injections for the first mirror
    PTT2 = logical([1 1 1]); % turns on/off different piston, tip/tilt difference field injections
    
    if looping == true
        %     [TipEstimate] = input('Enter appropriate tip amount: ');
        PTTpositions1 = [Estimate{1} 1e-3 1e-3]; % default piston tip/tilt positions for the difference field segment for the first dOTF image
        PTTpositions2 = [lambda/4 1e-3 1e-3]; % default piston tip/tilt positions for the difference field segment for the second dOTF image
    else
        PTTpositions1 = [lambda/8 1e-3 1e-3]; % default piston tip/tilt positions for the difference field segment for the first dOTF image
        PTTpositions2 = [lambda/4 1e-3 -1e-3]; % default piston tip/tilt positions for the difference field segment for the second dOTF image
    end
    
    for i = 1:length(PTT1)
        if isempty(PTTpositions1(PTT1(i)))
            PTTpositions1(i) = 0;
        end
        
        if isempty(PTTpositions2(PTT2(i)))
            PTTpositions2(i) = 0;
        end
    end
    
    if isequal(PTTpositions1,PTTpositions2)
        error('Your dOTF is nonexistent! Give an actual difference.');
    end
    
    PTTpos_poked1(DiffField,1:3) = PTTpositions1;
    PTTpos_poked2(DiffField,1:3) = PTTpositions2;
    
    % Map JLC's segment numbering to IrisAO segment numbering
    PTT_poked1 = mapSegments(PTTpos_poked1,numRings);
    PTT_poked2 = mapSegments(PTTpos_poked2,numRings);
    
else
    PTT_poked1 = mapSegments(PTTpos_poked1,numRings);
    PTT_poked2 = mapSegments(PTTpos_poked2,numRings);
end

% Generate the exact pupil modification
DMDiff.PTT(PTTpositions1);
DMDiff.touch;
DMDiff.render;

FD.planewave * ABER * TURB * DMDiff;
% FD.show;
fdiffexact = FD.grid;
FDIFFEXACT = fftshift(fft2(fftshift(fdiffexact)));

DMDiff.PTT([0 0 0]);
DMDiff.touch;
DMDiff.render;
FD.planewave * DMDiff;
% FD.show;

fdiffoffset = FD.grid;

fdiffmeasured = fdiffexact - fdiffoffset;
FDIFFMEASURED = fftshift(fft2(fftshift(fdiffmeasured)));

% Render the modified mirror surface
DM1.PTT(PTT_poked1);
DM1.touch;
DM1.render;
% figure;
% DM1.show; colorbar;

F.planewave * ABER * TURB * A * DM1;
% F.planewave * ABER * TURB * DM1;

if segment ~= true
    [XF, YF] = F2.COORDS;
    F2.grid(abs(YF-14.3e-4) < .0051/40 & abs(XF-8.7e-4) < (0.0051/80));
    % F2.grid(exp(1i*k*double(abs(YF-14.3e-4) < .0051/40 & abs(XF-8.7e-4) < (0.0051/80))));
    (F - F2);
end

f1 = F.grid; % save for computing the difference field


DM1.PTT(PTT_poked2);
DM1.touch;
DM1.render;
% figure;
% DM1.show; colorbar;

F2.planewave * ABER * TURB * A * DM1;
% F2.planewave * ABER * TURB * DM1;

f2 = F2.grid; % save the update to the mirror surface

%% dOTF

% Generate Difference field
if exactDifference == true
    fdiff = fdiffmeasured;
    FDIFF = fftshift(fft2(fftshift((fdiff))));
else
    fdiff = f1 - f2;
    load P.mat
    if manualPick == true
        figure;
        plotComplex(fdiff,2);
        P = pickPoint;
        if pointOverwrite == true
            save P.mat P
        end
    end
    FDIFF = fftshift(fft2(circshift(conj(fdiff),1-P)));
end
% FD.grid((circshift(FD.grid, 1 - P))); % place at corners

% set up comparison case
halo = fftshift(fft2(fftshift(F.grid)));
PupilP = fftshift(ifft2(ifftshift(halo.*FDIFF)));
%     figure;
%     plotComplex(PupilP,6); axis xy; sqar;


PSF1 = abs(fftshift(fft2(fftshift(f1)))).^2;
PSF2 = abs(fftshift(fft2(fftshift(f2)))).^2;

if Noise_Parameters{5} == true
    [PSF1,PSF2] = average_noisy_images(PSF1,PSF2,Noise_Parameters{6},Noise_Parameters);
end

OTF1 = fftshift(fft2(fftshift(PSF1)));
OTF2 = fftshift(fft2(fftshift(PSF2)));

% Do the dOTF
dOTF = OTF2 - OTF1;
mag = abs(dOTF);
phase = angle(dOTF);
% unwrapped_phase = uwrap(phase,'unwt');

%% Masking
% Generating coordinates for masking dOTF
G = AOGrid(1);
G.spacing(DM1.spacing);
G.grid(dOTF);
[X, Y] = G.COORDS;
GT = G.copy;

% based on difference field segment location from the center segment
if Minus == true
    if smoothing == false
        GT.grid((Y - tand(150)*X > 0)); % fftsize = 2^10, manual ffts
    elseif smoothing == true
        GT.grid(smoothedge(Y - tand(150)*X - .2e-3,10*SPACING));
    end
    load POMinus.mat
else
    if smoothing == false
        GT.grid(~(Y - tand(150)*X > 0));
    elseif smoothing == true
        GT.grid(smoothedge(-(Y - tand(150)*X + .2e-3),10*SPACING));
    end
    load PT.mat
end

% Apply the masking
G * GT;

if manualPick == true
    figure; plotComplex(G.grid,2); axis xy; sqar;
    PT = pickPoint;
    if pointOverwrite == true && Minus == true
        save POMinus.mat PT
    elseif pointOverwrite == true && Minus == false
        save PT.mat PT
    end
end

G.grid(fftshift(circshift(G.grid, 1 - PT))); % centers the dOTF grid

%     maskD = 3e-3; % 2 mm
%     PUPIL_DEFN = [
%         0 0 maskD         1 aa 0 0 0 0 0
%         %         0 0 secondary 0 aa/2 0 0 0 0 0
%         %         0 0 spider   -2 aa 4 0 D/1.9 0 0
%         ];
%     M = AOSegment;
%     M.spacing(SPACING);
%     M.name = 'dOTF Mask';
%     M.pupils = PUPIL_DEFN;
%     M.make;
%     mask.planewave * M;
%     G * mask;


%% Deconvolution
% Fourier Regularizated Deconvolution without noise statistics right now
DOTF = fftshift(fft2(fftshift(G.grid)));
TDOTF = fftshift(fft2(fftshift(dOTF)));
gamma = .7e4; % regularization parameter
if Minus == true
    Deconv = Wiener(DOTF,FDIFF,1,gamma);
else
    Deconv = Wiener(DOTF,conj(FDIFF),1,gamma);
end
DeconvP = Wiener(fftshift(fft2(fftshift(PupilP))),FDIFF,1,gamma);
DeconvTotal = Wiener(TDOTF,FDIFF,1,gamma);


%% Debugging
%Mostly just plotting figures right now to compare intermediate steps in
%the dOTF calculation.
if debuggery == true && Plotting == true
    
    M1 = 2; M2 = 4;
    fontsize = 10;
    
    figure;
    %     subplot(M1,M2,1);
    %     imagesc(log10([PSF1 PSF2]));
    %     axis xy; sqar;
    %     bigtitle('PSF1 Left, PSF2 Right',fontsize);
    
    subplot(M1,M2,1);
    plotComplex([OTF1 OTF2],2);
    axis xy; sqar;
    bigtitle('OTF1 Left, OTF2 Right',fontsize);
    
    subplot(M1,M2,2);
    plotComplex(fftshift(fft2(fftshift(DeconvTotal))),6);
    axis xy; sqar;
    bigtitle('Fourier tranform of the deblurred dOTF',fontsize);
    
    subplot(M1,M2,3);
    plotComplex(FDIFF,6);
    axis xy; sqar;
    bigtitle('Fourier transform of the difference field',fontsize);
    
    subplot(M1,M2,4);
    plotComplex(FDIFFEXACT,6);
    axis xy; sqar;
    bigtitle('Fourier transform of the exact difference field',fontsize);
    
    subplot(M1,M2,5);
    plotComplex(FDIFFMEASURED,6);
    axis xy; sqar;
    bigtitle('Difference Field Generation Test',fontsize);
    
    subplot(M1,M2,6);
    plotComplex(fftshift(fft2(fftshift(Deconv))),6);
    axis xy; sqar;
    bigtitle('Fourier transform of the Masked dOTF',fontsize);
    
    subplot(M1,M2,7);
    plotComplex(halo,7);
    axis xy; sqar;
    bigtitle('Fourier transform of the Pupil',fontsize);
    
    subplot(M1,M2,8);
    plotComplex(halo.*FDIFF,6);
    axis xy; sqar;
    bigtitle('Fourier transform of the Blurred Pupil',fontsize);
    
end
%% Plotting
if Plotting == true
    N1 = 2; N2 = 4; fontsize = 10;
    figure;
    
    subplot(N1,N2,1)
    plotComplex(dOTF,2);
    axis xy; sqar;
    bigtitle('Complex Plot of dOTF',fontsize);
    
    subplot(N1,N2,2)
    imagesc(mag);
    colormap(gray); axis xy; colorbar; sqar;
    bigtitle('Magnitude of dOTF',fontsize);
    
    subplot(N1,N2,3)
    imagesc(phase);
    colormap(gray); axis xy; colorbar; sqar;
    bigtitle('Phase of dOTF',fontsize);
    
    subplot(N1,N2,4)
    plotComplex(PupilP,2);
    axis xy; sqar;
    bigtitle('Blurred Pupil',fontsize);
    
    subplot(N1,N2,5)
    plotComplex(G.grid,2);
    axis xy; sqar;
    bigtitle('Masked and centered dOTF Pupil Estimate',fontsize);
    
    subplot(N1,N2,6)
    if Minus == true
        plotComplex(rot90(-Deconv,2),2);
    else
        plotComplex(-conj(Deconv),2);
    end
    axis xy; sqar;
    bigtitle('dOTF Deconvolution',fontsize);
    
    subplot(N1,N2,7)
    plotComplex(rot90(-DeconvTotal,2),2);
    axis xy; sqar;
    bigtitle('Total Deconvolved dOTF',fontsize);
    
    subplot(N1,N2,8)
    plotComplex(DeconvP,2);
    axis xy; sqar;
    bigtitle('Deblurred Pupil',fontsize);
    
    % imagesc(unwrapped_phase .* mask,[-2*pi,2*pi]);
    % colormap(gray); axis xy; colorbar; sqar;
end

if looping == true
    if Estimate{2} > Estimate{3}
        Estimate{1} = 1e-4;
        Estimate{2} = 1;
        return;
    else
        Estimate{1} = Estimate{1} + lambda/64;
        Estimate{4}(:,:,Estimate{2}) = DeconvTotal;
        Estimate{2} = Estimate{2} + 1;
        Estimate = dOTF_Sim_Builder(Estimate);
    end
end
% end
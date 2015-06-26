clear all;
clc;
close all;

% Closed-Loop dOTF Methods and Testing

%**************************************************************************
%                       SPIE dOTF Closed-Loop Methods Simulation
%**************************************************************************

%% System Parameters
% Don't touch these if using maketestbedelements.m to construct components
global lambda k D secondary spider SPACING aa fftsize THld FOV PLATE_SCALE FoV_withIrisAO FoV_withoutIrisAO RunSIM RunTESTBED IrisAO_on BMC_on verbose_makeDM Scalloped_Field UseRealPSF coronagraph system_verbose

% Set Wavelength
% lambda = AOField.RBAND; % Red light.
lambda = AOField.HeNe_Laser;

% Compute Wavenumber
k = (2*pi) / lambda;

% Pupil Specs
D = 3.636e-3; % 7mm
% secondary = 0.3*D; % 30% of Pupil Diameter
secondary = 0;
% spider = 0.02*D; % 2% of Pupil Diameter
spider = 0;


%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 4; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
goal_strehl = 0.9; %exit condition for loop
gain = 0.7; %gain for AO Corrections
fftsize = 2^12;

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 100*THld; % arcsecs
PLATE_SCALE = THld/4;
FoV_withIrisAO = 3.5e-3;
FoV_withoutIrisAO = 10.5e-3;

%% Run Simulation/Testbed Flags
RunSIM = true; %Run the simulation
RunTESTBED = false; %Run the testbed equipment

%% Simulation Flags

% IrisAO Flags
IrisAO_on = true; %turns on/off IrisAO Mirror (if false, DM1 variable set to 1)
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = false; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.
spiral = true;

% BMC Flag
BMC_on = false; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = false; %Injects nzerns Zernike Terms
InjectRandAb = false; %if InjectAB is true, picks Zernikes "Randomly"
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
    number_of_images = 100;
    Noise = cell(2,1);
    Noise{1} = UseNoise;
    Noise{2} = number_of_images;
else
    Noise = cell(1,1);
    Noise{1} = UseNoise;
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
    varargin{1} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithoutFingerDMBox/';
    varargin{3} = 'RAW_scienceIM_frame_';
    varargin{2} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithFingerDMBox/';
    varargin{4} = 'RAW_scienceIM_frame_';
end
% Coronagraph Flag
coronagraph = false; % turns on going through coronagraph elements

% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Flags



%% Make the Testbed Elements
numRings = 3;
numSeg = sum(1:numRings)*6 + 1;
load MadeTestbedElements.mat;
% maketestbedelements;


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
        
        TURB = AOAtmo(A);
        %         TURB.spacing(SPACING);
        WFlow = AOScreen(fftsize,0.15,500e-9);
        WFlow.name = 'Lower altitude turbulence';
        WFhigh = AOScreen(2*fftsize,0.17,500e-9);
        WFhigh.name = 'High altitude turbulence';
        
        TURB.addLayer(WFlow,1000);
        TURB.addLayer(WFhigh,8000);
        
        TURB.layers{1}.Wind = [3 1];
        TURB.layers{2}.Wind = [1 -1]*20;
        
        r0 = TURB.totalFriedScale;
        th_scat = lambda/r0*206265;
        
        fprintf('The total r0 is %f cm.\n',100*r0);
        fprintf('The seeing is %.2f arcsecs.\n',th_scat);
        
        % Turning this off is like using dynamic refocus.
        TURB.GEOMETRY = true;
        TURB.BEACON = [1 1 1e10];
        TURB.make;
        TURB.show
        
        aberration = TURB.grid;
        
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

% Create the location of the difference field
if numRings > 0
    DiffField = 3;
    for segRing = 1:numRings-1
        DiffField = DiffField + (segRing)*6 + 1;
    end
else
    DiffField = 1;
end

% Set Finger Positions
PTTpos_poked2 = PTTpos_poked1;
PTTpos_poked1(DiffField,1) = 1e-6; % 1 micron
% PTTpos_poked1(DiffField,2) = 1e-4; % .1 mrad
PTTpos_poked2(DiffField,3) = 1e-4; % .1 mrad
PTT_poked1 = mapSegments(PTTpos_poked1,numRings);
PTT_poked2 = mapSegments(PTTpos_poked2,numRings);

DM1.PTT(PTT_poked1);
DM1.touch;
DM1.render;
% figure;
% DM1.show; colorbar;

F.planewave * ABER * TURB * A * DM1;
% figure; F.show;
PSF1 = F.mkPSF(FOV,PLATE_SCALE);
PSF1max = max(max(PSF1));
% PSF1 = PSF1 / PSF1max;
PSF1plot = log10(PSF1/PSF1max);
% figure; imagesc(PSF1plot);

DM1.PTT(PTT_poked2);
DM1.touch;
DM1.render;
% figure;
% DM1.show; colorbar;

F2.planewave * ABER * TURB * A * DM1;
PSF2 = F2.mkPSF(FOV,PLATE_SCALE);
PSF2max = max(max(PSF2));
% PSF2 = PSF2 / PSF2max;
PSF2plot = log10(PSF2/PSF2max);

% figure(1);
% imagesc([PSF1plot,PSF2plot],[-4,0]);
% colormap(gray);

F.touch; F2.touch;
F.grid(PSF1); F2.grid(PSF2);

OTF1 = F.mkOTF2(FoV_withoutIrisAO,SPACING(1));
OTF2 = F2.mkOTF2(FoV_withoutIrisAO,SPACING(1));

F.touch; F2.touch;

% figure(2)
% subplot(1,2,1)
% plotComplex(OTF1,2);
% subplot(1,2,2)
% plotComplex(OTF2,2);

dOTF = OTF1 - OTF2;
mag = abs(dOTF);
phase = angle(dOTF);
% unwrapped_phase = uwrap(phase,'unwt');

% mask = mag;
% mask(mask<1e10) = 0;
% mask = logical(mask);


figure(3);
subplot(1,3,1)
% imagesc(mag .* mask);
imagesc(mag);
colormap(gray); axis xy; colorbar; sqar;
subplot(1,3,2)
% imagesc(phase .* mask);
imagesc(phase);
colormap(gray); axis xy; colorbar; sqar;
subplot(1,3,3)
plotComplex(dOTF,2);
axis xy; sqar;
% imagesc(unwrapped_phase .* mask,[-2*pi,2*pi]);
% colormap(gray); axis xy; colorbar; sqar;

















































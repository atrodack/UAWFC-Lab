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
D = 7e-3; % 7mm
secondary = 0.3*D; % 30% of Pupil Diameter
spider = 0.02*D; % 2% of Pupil Diameter


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
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.

% BMC Flag
BMC_on = true; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = true; %Injects nzerns Zernike Terms
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
        if strcmpi(answer1,'y') == false
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
coronagraph = false; % turns on going through coronagraph elemens

% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Flags



%% Make the Testbed Elements
maketestbedelements;
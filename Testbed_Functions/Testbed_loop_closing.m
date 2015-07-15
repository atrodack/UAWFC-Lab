clear all;
clc;
% close all;

% Closed-Loop dOTF Methods and Testing

%**************************************************************************
%                       SPIE dOTF Closed-Loop Methods Simulation
%**************************************************************************

%% System Parameters

% Set Wavelength
% lambda = AOField.RBAND; % Red light.
lambda = AOField.HeNe_Laser;

% Compute Wavenumber
k = (2*pi) / lambda;

% Pupil Specs
D = 7e-3; % 7mm
% secondary  0.3*D; % 30% of Pupil Diameter
secondary = 0;
% spider = 0.02*D; % 2% of Pupil Diameter
spider = 0;


%% Simulation Flags
% Use Testbed PSF Instead of Simulated PSF
UseRealPSF = true;
if UseRealPSF == true
    % Num_Folders == Number of folders images are to be taken from
    Num_Folders = 2;
    
    % Num_files_per_folder == Number of images in the folder to load in
    Num_files_per_folder = 100;
    
    % varargin{1-Num_Folders} == path to folder 1, varargin{Num_Folders+1 -2*Num_Folders} == name of individual image file sans number at end
    varargin{1} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithoutFingerDMBox/';
    varargin{3} = 'RAW_scienceIM_frame_';
    varargin{2} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithFingerDMBox/';
    varargin{4} = 'RAW_scienceIM_frame_';
end

% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Stuff
pokeact = 698;



%% Make DM Model
makeBMCmirror;




%% Write a Filename
dt = datestr(now,'mm_dd_yyyy_HH_MM_SS_FFF');
filename = sprintf('Testbed_BMC_Closed_loop_%s',dt);
filename_movie_avi = sprintf('%s.avi',filename);
fprintf('\n');
display(filename);
fprintf('\n');



%% Calibrate the Pixel-Actuator Map (if not done/saved somewhere)
% calibrated_BMC_act_locations = computetestbedBMCactpixelmap(DM,pokeact);













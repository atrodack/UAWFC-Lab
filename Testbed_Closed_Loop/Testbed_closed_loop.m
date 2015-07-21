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


%% Testbed Stuff
DISK = fitsread('DISK.fits');
DM_Acts = DISK(:);

actlist = zeros(length(DM_Acts),2);
actlist(:,1) = linspace(1,length(DM_Acts),length(DM_Acts));
actlist(:,2) = DM_Acts;

% Create a mapping to the registered actuators in DISK
% Coordinates in DISK
[Act_map(:,1),Act_map(:,2)] = find(abs(DISK)>0);
% Corresponding Actuator Number
Act_map(:,3) = actlist(actlist(:,2)~=0);

% The ordering is (as noted in the first two columns) actuator number 1 is
% bottom left, and increases by ascending the column, moving down to up,
% left to right

Off_Acts = actlist(actlist(:,2) ==0);



Run_Cam_Parameters{1} = 100;
Run_Cam_Parameters{2} = 48;
Run_Cam_Parameters{3} = false;
Run_Cam_Parameters{4} = 1;
Run_Cam_Parameters{5} = 'HeNe';
Run_Cam_Parameters{6} = 'OD3';


%% Make DM Model
load BMC_DM_Model.mat

% Calibrate the Model
DM.disableActuators(Off_Acts);
DM.setOnActs;
DM.enableActuators(Off_Acts);
DM.settheOnActs(DM_Acts(Act_map(:,3)));

% The Model is set to what one would see in the shmimview windows
% ie not with the 90 degree rotation

%% Write a Filename
dt = datestr(now,'mm_dd_yyyy_HH_MM');
filename = sprintf('Testbed_Data_%s',dt);
filename_movie_avi = sprintf('%s.avi',filename);
fprintf('\n');
display(filename);
fprintf('\n');



%% Calibrate the Pixel-Actuator Map (if not done/saved somewhere)
% calibrated_BMC_act_locations = computetestbedBMCactpixelmap(DM,pokeact);
Ppos_in = fitsread('Testpattern1.fits');


% Do a dOTF with the testpattern
[dOTF_cal, PSF_CUBE_cal, PSF_poked_CUBE_cal] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in);

mkdir(filename);
current_dir = pwd;
cd(filename)
save('Data_Cubes_cal.mat','PSF_CUBE_cal','PSF_poked_CUBE_cal');
cd(current_dir);

phase = angle(dOTF_cal);
uphase = uwrap(phase,'unwt');
OPL = uphase / k;



%% Plot
figure;
subplot(1,3,1);
imagesc(PSF_CUBE_cal.PSF_centered_and_cropped); axis xy; axis off; sqar; bigtitle('PSF',15); %colormap(gray);

subplot(1,3,2)
plotComplex(dOTF_cal,5); axis xy; axis off; sqar; bigtitle('dOTF',15);

subplot(1,3,3);
imagesc(OPL); axis xy; axis off; sqar; bigtitle('OPL',15);



%%
% [X,Y] = meshgrid(1:1:256);
% R = sqrt((X-128).^2 + (Y-152).^2);
% mask = double(R<=26);
% 
% masked_OPL = mask.*OPL;
% masked_OPL = fftshift(circshift(masked_OPL, 1- [152,128]));
% OPL_ = masked_OPL(129-27:129+27,129-27:129+27);
% OPL_ = rot90(OPL_,1);
% OPL_binned = downsampleCCD(OPL_,3,3);
% DM_shape = (padarray(OPL_binned,[(32-18)/2,(32-18)/2]));






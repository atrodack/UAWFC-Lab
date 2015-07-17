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
pokeact = 698;



%% Make DM Model
load BMC_DM_Model.mat




%% Write a Filename
dt = datestr(now,'mm_dd_yyyy_HH_MM_SS_FFF');
filename = sprintf('Testbed_BMC_Closed_loop_%s',dt);
filename_movie_avi = sprintf('%s.avi',filename);
fprintf('\n');
display(filename);
fprintf('\n');



%% Calibrate the Pixel-Actuator Map (if not done/saved somewhere)
% calibrated_BMC_act_locations = computetestbedBMCactpixelmap(DM,pokeact);


%% Do a dOTF
[dOTF, PSF_CUBE, PSF_poked_CUBE] = TestbeddOTF(DM, pokeact, false);



phase = angle(dOTF);
uphase = uwrap(phase,'unwt');
OPL = uphase / k;



%% Plot
figure(1);
subplot(1,3,1);
imagesc(PSF_CUBE.PSF_dark_corr); axis xy; axis off; sqar; bigtitle('PSF',15); %colormap(gray);

subplot(1,3,2)
plotComplex(dOTF,5); axis xy; axis off; sqar; bigtitle('dOTF',15);

subplot(1,3,3);
imagesc(OPL); axis xy; axis off; sqar; bigtitle('OPL',15);





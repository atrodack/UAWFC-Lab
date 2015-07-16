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


%% Do a dOTF
% [dOTF, PSF, PSF_poked, OTF, OTF_poked] = TestbeddOTF(DM, pokeact, false);
[PSF_CUBE, PSF_poked_CUBE] = TestbeddOTF(DM, pokeact, false);



% phase = angle(dOTF);
% uphase = uwrap(phase,'unwt');
% OPL = uphase / k;

%% Plot the Results of the dOTF
% figure(1);
% subplot(2,3,1)
% imagesc(PSF); axis xy; sqar; colormap(gray); bigtitle('PSF',15);
% 
% subplot(2,3,4);
% imagesc(PSF_poked); axis xy; sqar; bigtitle('PSF with Pupil Mod',15);
% 
% subplot(2,3,2);
% plotComplex(OTF,6); axis xy; sqar; bigtitle('OTF',15);
% 
% subplot(2,3,5);
% plotComplex(OTF_poked,6); axis xy; sqar; bigtitle('OTF with Pupil Mod',15);
% 
% subplot(2,3,3)
% plotComplex(dOTF,6); axis xy; sqar; bigtitle('dOTF',15);
% 
% subplot(2,3,6)
% imagesc(OPL); axis xy; sqar; axis off; colorbar; bigtitle('OPL',15);








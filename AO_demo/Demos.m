clear all;
clc;
closepreview;
% close all;
% cd ~/Desktop/UAWFC/ % home directory

%==========================================================================
%                         Closed-Loop dOTF Demos
%==========================================================================

% folders on MATLAB pathway ~/Desktop/UAWFC/: AO_demo, AOSim2, dOTF_Calibrations, Testbed_Functions


%% System Parameters
% Pupil Specs
D = 7e-3; % 7mm
secondary = 0.3*D; % 30% of Pupil Diameter
spider = 0.02*D; % 2% of Pupil Diameter

%% Testbed Parameters
% Figure out the illuminated actuators on the DM using previously
% calibrated locations stored in fits file.

% Figure out illuminated actuators using calibrated beam position
DISK = fitsread('DISK_DOWN.fits');
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

% Make a list of the actuators that are not illuminated
Off_Acts = actlist(actlist(:,2) ==0);

% Initial Mirror Position
Ppos_in = zeros(32,32);

% Command to zero poke channel
cmd_clear_poke = sprintf('~/src/scripts/dmzeroch %d',0);

%% Make DM Model and set actuators to above list
% The functions run default to using actuator lists from this model if no
% inputs are supplied. It might not be necessary to include, but for now it
% does...
% The Model is set to what one would see in the shmimview windows
% i.e. not with the 90 degree rotation
load BMC_DM_Model.mat

% Calibrate the Model
DM.disableActuators(Off_Acts);
DM.setOnActs;
DM.enableActuators(Off_Acts);
DM.settheOnActs(DM_Acts(Act_map(:,3)));




%% Write a Filename using the Date
dt = datestr(now,'mm_dd_yyyy_HH_MM');
filename = sprintf('Testbed_Data_%s',dt);
% filename_movie_avi = sprintf('%s.avi',filename);
fprintf('\n');
display(filename);
fprintf('\n');



%% Calibrate the Pixel-Actuator Map (if not done/saved somewhere)
% Must calibrate the location of the DM actuators with the dOTF pupil
% measurement.
% Already done for using actuator [17,8] as poke actuator, make sure this
% is set for doing closed-loop


%% dOTF Masks (this might be deprecated...)
% Make masks for obtaining the correct OPL from dOTF measurements
% Also made for using actuator [17,8]
[X,Y] = meshgrid(1:1:256);
R = sqrt((X-129).^2+ (Y-103).^2);
mask1 = double(R<=30);
R2 = sqrt((X-129).^2+ (Y-155).^2);
mask2 = double(R2<= 30);
mask = mask1 + mask2;
mask(mask~=0) = 1;

%% Demonstration #1
% Exploring ways to obtain dOTF measurements
%--------------------------------------------------------------------------

%% 1.	Amplitude Modification (CHANGE CENTERPOINT)
% If a DM is not an element in the optical system for which a phase
% measurement is desired, a simple amplitude blocker can be used. This can
% be anything from a thin wire to a small piece of tape. Images are taken
% at the science camera with no blocker to obtain the unmodified pupil OTF.
% Then, the blocker is placed near the edge of the illuminated pupil within
% a pupil plane of the optical system, and images are taken again. Where the
% blocker is placed in the pupil determines the geometry of the conjugates in
% the dOTF measurement.
%
% Observe that as the blocker is moved in the pupil, the orientation of the
% pupils in the dOTF changes, as well as the size of the overlap region.
%==========================================================================

% Set Wavelength
% lambda = 600*10^-9; % Central wavelength of bandpassed SuperK
lambda = AOField.HeNe_Laser;

% Band Pass Filter
bandpassfilter = 0; %0=no filter, 40=40nm band, 10=10nm band

% Compute Wavenumber
k = (2*pi) / lambda;

% Set Camera Parameters
Run_Cam_Parameters{1} = 5; %Number of Triggers (0 based, want 10 triggers, this equals 9)
Run_Cam_Parameters{2} = 30; %Number of Frames per Trigger
Run_Cam_Parameters{3} = 50; %Shutter Speed of Camera
Run_Cam_Parameters{4} = false; %Take Darks: T/F
Run_Cam_Parameters{5} = 1; %Camera Channel
Run_Cam_Parameters{8} = 16; %Gain
Run_Cam_Parameters{9} = 60; %Exposure

if lambda == AOField.HeNe_Laser
    Run_Cam_Parameters{6} = 'HeNe'; %Laser Type
    Run_Cam_Parameters{7} = 'None'; %Filter Type
else
    Run_Cam_Parameters{6} = 'SuperK';
    if bandpassfilter == 0
        Run_Cam_Parameters{7} = 'None'; %Filter Type
    elseif bandpassfilter == 40
        Run_Cam_Parameters{7} = '40 nm'; %Filter Type
    elseif bandpassfilter == 10
        Run_Cam_Parameters{7} = '10 nm'; %Filter Type
    end
end


% Take images with no blocker in place
fprintf('\nTaking images with no blocker\n');
PSF_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);
fprintf('\nNo blocker PSFs Complete\n');

% Remind to put in blocker
input('\nPress Enter when blocker is in place\n');

% Take images with the blocker in
PSF_poked_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);
fprintf('\nBlocker PSFs Complete\n');

% Process
[PSF_CUBE] = ExtractAverageCUBEPSFs_AO_demo(PSF_CUBE);
[PSF_poked_CUBE] = ExtractAverageCUBEPSFs_AO_demo(PSF_poked_CUBE);
[dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF2_AO_demo(PSF_CUBE, PSF_poked_CUBE);
dOTF(129,:) = 0;
dOTF = -1i * conj(dOTF);

figure(2);
subplot(1,3,1)
imagesc(PSF_CUBE.PSF_centered_and_cropped);
axis xy;
sqar;
title('PSF without Blocker');

subplot(1,3,2)
imagesc(PSF_poked_CUBE.PSF_centered_and_cropped);
axis xy;
sqar;
title('PSF with Blocker');

subplot(1,3,3)
plotComplex(dOTF,2);
axis xy;
sqar;
title('dOTF');
drawnow;

%% 2.	Phase Modification (CHANGE CENTERPOINT)
% If a DM is an element in the optical system and is located in a pupil plane,
% an actuator can be used as the pupil modification. Using this method over
% an amplitude blocker gives us another degree of freedom: the magnitude of
% the actuator poke. We will now explore what happens to the signal for a
% couple of different magnitudes of phase modifications.
%
% Notice that there will be an optimal magnitude poke to use.

% Set Wavelength
lambda = 600*10^-9; % Central wavelength of bandpassed SuperK
% lambda = AOField.HeNe_Laser;

% Band Pass Filter
bandpassfilter = 0; %0=no filter, 40=40nm band, 10=10nm band

% Compute Wavenumber
k = (2*pi) / lambda;

% Set Camera Parameters
Run_Cam_Parameters{1} = 4; %Number of Triggers (0 based, want 10 triggers, this equals 9)
Run_Cam_Parameters{2} = 30; %Number of Frames per Trigger
Run_Cam_Parameters{3} = 40; %Shutter Speed of Camera
Run_Cam_Parameters{4} = false; %Take Darks: T/F
Run_Cam_Parameters{5} = 1; %Camera Channel
Run_Cam_Parameters{8} = 16; %Gain
Run_Cam_Parameters{9} = 60; %Exposure

if lambda == AOField.HeNe_Laser
    Run_Cam_Parameters{6} = 'HeNe'; %Laser Type
    Run_Cam_Parameters{7} = 'None'; %Filter Type
else
    Run_Cam_Parameters{6} = 'SuperK';
    if bandpassfilter == 0
        Run_Cam_Parameters{7} = 'None'; %Filter Type
    elseif bandpassfilter == 40
        Run_Cam_Parameters{7} = '40 nm'; %Filter Type
    elseif bandpassfilter == 10
        Run_Cam_Parameters{7} = '10 nm'; %Filter Type
    end
end

% Set magnitude of actuator pokes
pokeamps = [(lambda*10^6), (lambda*10^6)/2, (lambda*10^6)/4];


for nn = 1:length(pokeamps)
    
    [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],pokeamps(nn));
    dOTF(129,:) = 0;
    dOTF = -1i * conj(dOTF);
    
    figure(2);
    subplot(1,3,nn)
    plotComplex(dOTF,2);
    axis xy;
    sqar;
    title(sprintf('dOTF with\n%0.3f micron poke\n',pokeamps(nn)));
    drawnow;
    
    if nn < length(pokeamps)
%         input('\nPress Enter to move on to next actuator magnitude\n');
    end
end %for
fprintf('\n\n');
system(cmd_clear_poke);
%% 3.	Wavelength Dependence
% The dOTF depends on the wavelength of the source being used. With all
% other things being equal in each dOTF measurement, the narrower the
% bandwidth of the source, the less blurred the measurement will be.
% We can observe this by changing the bandpass filter of the white light laser source.


% Band Pass Filter
bandpassfilter = [10, 40, 0]; %0=no filter, 40=40nm band, 10=10nm band
sSpeed = [40,3,1];

% NOTE:
% OD filter must be removed for 10nm bandpass

for nn = 1:3
    % Set Camera Parameters
    Run_Cam_Parameters{1} = 5; %Number of Triggers (0 based, want 10 triggers, this equals 9)
    Run_Cam_Parameters{2} = 30; %Number of Frames per Trigger
    Run_Cam_Parameters{3} = sSpeed(nn); %Shutter Speed of Camera
    Run_Cam_Parameters{4} = false; %Take Darks: T/F
    Run_Cam_Parameters{5} = 1; %Camera Channel
    Run_Cam_Parameters{8} = 16; %Gain
    Run_Cam_Parameters{9} = 60; %Exposure
    
    if lambda == AOField.HeNe_Laser
        Run_Cam_Parameters{6} = 'HeNe'; %Laser Type
        Run_Cam_Parameters{7} = 'None'; %Filter Type
    else
        Run_Cam_Parameters{6} = 'SuperK';
        if bandpassfilter(nn) == 0
            Run_Cam_Parameters{7} = 'None'; %Filter Type
        elseif bandpassfilter(nn) == 40
            Run_Cam_Parameters{7} = '40 nm'; %Filter Type
        elseif bandpassfilter(nn) == 10
            Run_Cam_Parameters{7} = '10 nm'; %Filter Type
        end
    end
    
    
    [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
    dOTF(129,:) = 0;
    dOTF = -1i * conj(dOTF);
    
    figure(4);   
    subplot(1,3,nn)
    plotComplex(dOTF,3);
    axis xy;
    sqar;
    if nn ~= 3
        title(sprintf('dOTF with\n%d nm bandpass filter\n',bandpassfilter(nn)));
    else
        title(sprintf('dOTF with\n no bandpass filter\n'));
    end
    drawnow;
    
    if nn < 3
        input('\nPress Enter after filter is changed\n');
    end
end %for
system(cmd_clear_poke);

%% Demonstration #2
% Now we will look more carefully at how the dOTF technique can be used as
% a wavefront sensor. In the previous demonstrations, it should have been
% obvious that there is some residual aberration in this system because the
% dOTF measurements were not constant across the pupil. This residual is
% due to non-common path aberrations and very slight misalignments. To make
% the power of the technique more clear, we will now add some low-order
% aberrations to the system using the DM and use the HeNe as the source.
%
% Try to guess what aberrations have added to the system by looking at the dOTF measurement that is displayed.

% Set Wavelength
lambda = 600*10^-9; % Central wavelength of bandpassed SuperK

% Band Pass Filter
bandpassfilter = 10; %0=no filter, 40=40nm band, 10=10nm band

% Compute Wavenumber
k = (2*pi) / lambda;

% Set Camera Parameters
Run_Cam_Parameters{1} = 5; %Number of Triggers (0 based, want 10 triggers, this equals 9)
Run_Cam_Parameters{2} = 30; %Number of Frames per Trigger
Run_Cam_Parameters{3} = 50; %Shutter Speed of Camera
Run_Cam_Parameters{4} = false; %Take Darks: T/F
Run_Cam_Parameters{5} = 1; %Camera Channel
Run_Cam_Parameters{8} = 16; %Gain - default 16
Run_Cam_Parameters{9} = 60; %Exposure

if lambda == AOField.HeNe_Laser
    Run_Cam_Parameters{6} = 'HeNe'; %Laser Type
    Run_Cam_Parameters{7} = 'None'; %Filter Type
else
    Run_Cam_Parameters{6} = 'SuperK';
    if bandpassfilter == 0
        Run_Cam_Parameters{7} = 'None'; %Filter Type
    elseif bandpassfilter == 40
        Run_Cam_Parameters{7} = '40 nm'; %Filter Type
    elseif bandpassfilter == 10
        Run_Cam_Parameters{7} = '10 nm'; %Filter Type
    end
end

% Zernike Low-order Aberration Files (fits files may or may not yet exist)
tempdir = pwd;
cd /home/lab/src/scripts/dOTF_Zerns
% Read in aberrations in case we want to plot them
Tilt = fitsread('Tip1.fits');
Defocus = fitsread('Defocus1.fits');
Astig = fitsread('Astig901.fits');
Coma = fitsread('Comax1.fits');
Spherical = fitsread('Spherical1.fits');
% Set DM Channels (phase mod channel is 0, DM correction from 4D is 7)
% Make all the same to overwrite, make different to combine
TiltChannel = 6;
DefocusChannel = 5;
AstigChannel = 4;
ComaChannel = 3;
SphericalChannel = 2;
% Write system command strings
cmd_Tilt = sprintf('~/src/scripts/dOTF_Zerns/dmloadch Tip1.fits %d',TiltChannel);
cmd_Defocus = sprintf('~/src/scripts/dOTF_Zerns/dmloadch Defocus1.fits %d',DefocusChannel);
cmd_Astig = sprintf('~/src/scripts/dOTF_Zerns/dmloadch Astig901.fits %d',AstigChannel);
cmd_Coma = sprintf('~/src/scripts/dOTF_Zerns/dmloadch Comax1.fits %d',ComaChannel);
cmd_Spherical = sprintf('~/src/scripts/dOTF_Zerns/dmloadch Spherical1.fits %d',SphericalChannel);
cd(tempdir);

%% SEND ABERRATIONS TO MIRROR
cd /home/lab/src/scripts/dOTF_Zerns
fprintf('\nSending Commands to DM\n\n');

% uncomment to send
% system(cmd_Tilt);
% system(cmd_Defocus);
system(cmd_Astig);
% system(cmd_Coma);
% system(cmd_Spherical);
cd(tempdir);

[dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
dOTF(129,:) = 0;
% dOTF(130,:) = 0;
% dOTF(:,129) = 0;
dOTF = -1i * conj(dOTF);

figure(5);
plotComplex(dOTF,2);
axis xy;
title('dOTF with Low-Order Aberration Added')
system(cmd_clear_poke);

%% Clear Aberration Channels
cmd_clear_TiltChannel = sprintf('~/src/scripts/dmzeroch %d',TiltChannel);
cmd_clear_DefocusChannel = sprintf('~/src/scripts/dmzeroch %d',DefocusChannel);
cmd_clear_AstigChannel = sprintf('~/src/scripts/dmzeroch %d',AstigChannel);
cmd_clear_ComaChannel = sprintf('~/src/scripts/dmzeroch %d',ComaChannel);
cmd_clear_SphericalChannel = sprintf('~/src/scripts/dmzeroch %d',SphericalChannel);

system(cmd_clear_TiltChannel);
system(cmd_clear_DefocusChannel);
system(cmd_clear_AstigChannel);
system(cmd_clear_ComaChannel);
system(cmd_clear_SphericalChannel);

%% Demonstration #3
% We will now show that using the dOTF technique enables an adaptive optics 
% system to sense its inherent aberrations and correct itself. Because the 
% technique uses images taken at the science camera, the encoded phase
% information includes any non-common path errors and mis-alignments through
% the whole system. By registering the location of the DM actuators to their
% corresponding locations in the dOTF measurement, wavefront piston and
% wavefront slopes can be read out to DM commands. Running a control loop
% using these DM updates allows the adaptive optics system to self-calibrate.
% 
% The following demonstration will correct the system aberration caused by
% the HeNe laser not hitting the first off-axis parabolic mirror (OAP) in 
% the optical system.

Testbed_closed_loop_AO_demo;



















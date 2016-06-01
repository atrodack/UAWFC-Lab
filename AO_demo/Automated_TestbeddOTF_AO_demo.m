function [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,Channel,poke_act,pokeamp)
% [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in,Channel, poke_act,pokeamp)
% Function to compute dOTF using the UAWFC testbed to obtain PSFs
%
% INPUTS:
% DM is an AODM model of a 1024 pin BMC Mirror
% Run_Cam_Parameters: See Automated_Testbed_run_cam_ver2.m
% Ppos_in is the initial actuator piston matrix
% Channel is the Cfits Channel to apply the update to
% poke_act is the pixel location in a 32x32 array of the actuator to poke
% poke_amp is the magnitude of the poke in microns
%
%
% OUTPUTS:
% dOTF is the differential Optical Transfer Function computed using the testbed PSFs
% PSF, PSF_poked are the averaged together, centered, and cropped PSFs


% NOTE:
% Check pauses in all scripts. Might be able to pick up a little time by
% shorting pauses

%% Input Argument Check
% Checking the Input Arguments
tic
if nargin < 2
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
    Ppos_in = reshape(Ppos_in,32,32);
    Channel = 0;
    poke_act = [17,8]; %Default, calibrated actuator
    pokeamp = (AOField.HeNe_Laser*10^6) / 4; %Default, lambda/4 for HeNe
    Run_Cam_Parameters{1} = 9;
    Run_Cam_Parameters{2} = 25;
    Run_Cam_Parameters{3} = 10;
    Run_Cam_Parameters{4} = false;
    Run_Cam_Parameters{5} = 1;
    Run_Cam_Parameters{6} = 'HeNe';
    Run_Cam_Parameters{7} = 'OD3';
    Run_Cam_Parameters{8} = 17;
    Run_Cam_Parameters{9} = 60;
    
elseif nargin < 3
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
    Channel = 0;
    poke_act = [17,8]; %Default, calibrated actuator
    pokeamp = (AOField.HeNe_Laser*10^6) / 4; %Default, lambda/4 for HeNe
elseif nargin < 5
    poke_act = [17,8]; %Default, calibrated actuator
    pokeamp = (AOField.HeNe_Laser*10^6) / 4; %Default, lambda/4 for HeNe
elseif nargin < 6
    pokeamp = (AOField.HeNe_Laser*10^6) / 4; %Default, lambda/4 for HeNe
elseif nargin > 6
    error('Number of Inputs is Wrong: dOTF = TestbeddOTF(DM, Run_Cam_Parameters, Ppos_in, Channel, poke_act)');
end

if ~isa(DM,'AODM')
    error('DM MUST BE an AODM Model of a 1024 Pin BMC Mirror: Class Type is Incorrect');
else
    if DM.nActs ~= 1024
        error('DM MUST BE an AODM Model of a 1024 Pin BMC Mirror: Number of Actuators in Model is Incorrect');
    end
end


if length(Ppos_in(:)) ~= 1024
    error('Ppos_in MUST HAVE 1024 Elements');
end


%% Parameter Initialization
% Initializing a Couple Variables
Ppos_flat = Ppos_in;


%% Initial Setup
% Remove corner pixels (not really actuators)
Ppos_flat(1,1) = 0; Ppos_flat(32,32) = 0; Ppos_flat(32,1) = 0; Ppos_flat(1,32) = 0;

% Set DM
Ppos_flat = single(Ppos_flat);

% Take the First PSF Images
tempdir = pwd;
cd /home/lab/src/scripts/
fitswrite(Ppos_flat,'DM_flat.fits');

cmd_flat = sprintf('~/src/scripts/dmloadch DM_flat.fits %d',Channel);
%SEND TO MIRROR
% input('Press Enter to Send to Mirror');
fprintf('\nSending Commands to DM\n\n');
system(cmd_flat);
cd(tempdir);


%CHECK FOR MIRROR UPDATE
% input('Press Enter once Mirror has Updated to Start Taking Images');
pause(2); % Probably can be reduced....

%TAKE IMAGES
PSF_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);


% input('Press Enter to Modify Pupil');
fprintf('\nModifying the Pupil\n\n');
% pause(1);
% ! ~/src/scripts/dmzeroch 0

%% Modified Setup
% Set DM and Take the Second PSF Image
Ppos_poked = Ppos_flat;

% Poke finger
Ppos_poked(poke_act(1),poke_act(2)) = pokeamp;

cd /home/lab/src/scripts
fitswrite(Ppos_poked,'DM_poked.fits');
cmd_poke = sprintf('~/src/scripts/dmloadch DM_poked.fits %d',Channel);

%SEND TO MIRROR
% input('Press Enter to Send to Mirror');
fprintf('\nSending Commands to the DM\n');
system(cmd_poke);
cd(tempdir);

%CHECK FOR MIRROR UPDATE
% input('Press Enter once Mirror has Updated to Start Taking Images');
pause(2); % Probably can be reduced.....

%TAKE IMAGES
PSF_poked_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);

% system(cmd_flat);

%% Process the Pictures

% USING NEW WAY IN ExtractdOTF2_AO_demo SHOULDN'T REQUIRE THESE STEPS
% [PSF_CUBE] = ExtractAverageCUBEPSFs_AO_demo(PSF_CUBE);
% [PSF_poked_CUBE] = ExtractAverageCUBEPSFs_AO_demo(PSF_poked_CUBE);

% Average Together Taken Images, Center/Crop PSFs, Compute OTFs and dOTF
% [dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF(PSF_CUBE, PSF_poked_CUBE,PT{1},PT{2});
% [dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF2(PSF_CUBE, PSF_poked_CUBE);

% Use the NEW WAY
[dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF2_AO_demo(PSF_CUBE, PSF_poked_CUBE,2);

toc

end
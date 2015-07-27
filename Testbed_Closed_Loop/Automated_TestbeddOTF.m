function [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in,Channel,PT)
% [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in,Channel)
% Function to compute dOTF using the UAWFC testbed to obtain PSFs
% 
% INPUTS:
% DM is an AODM model of a 1024 pin BMC Mirror
% Run_Cam_Parameters: See Automated_Testbed_run_cam_ver2.m
% Ppos_in is the initial actuator piston matrix
% Channel is the Cfits Channel to apply the update to
%
% OUTPUTS:
% dOTF is the differential Optical Transfer Function computed using the testbed PSFs
% PSF, PSF_poked are the averaged together, centered, and cropped PSFs
% OTF, OTF_poked are the OTFs computed from PSF and PSF_poked respectively



%% Input Argument Check
% Checking the Input Arguments

if nargin < 2
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
    Ppos_in = reshape(Ppos_in,32,32);
    Channel = 0;
    PT{1} = [];
    PT{2} = [];
Run_Cam_Parameters{1} = 9;
Run_Cam_Parameters{2} = 25;
Run_Cam_Parameters{3} = 100;
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
    PT{1} = [];
    PT{2} = [];
elseif nargin < 5
    PT{1} = [];
    PT{2} = [];
elseif nargin > 5
    error('Number of Inputs is Wrong: dOTF = TestbeddOTF(DM, Run_Cam_Parameters, Ppos_in, Channel)');
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
% Set DM and Take the First PSF Image
Ppos_flat(1,1) = 0; Ppos_flat(32,32) = 0; Ppos_flat(32,1) = 0; Ppos_flat(1,32) = 0;
Ppos_flat = single(Ppos_flat);
% 
% 
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
pause(2);
PSF_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);


% input('Press Enter to Modify Pupil');
fprintf('\nModifying the Pupil\n\n');
% pause(1);
% ! ~/src/scripts/dmzeroch 0

%% Modified Setup
% Set DM and Take the Second PSF Image
Ppos_poked = Ppos_flat;
Ppos_poked(17,7) = Ppos_poked(17,7) + (AOField.HeNe_Laser*10^6) / 4;


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
pause(2);
PSF_poked_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);

% system(cmd_flat);

%% Process the Pictures
[PSF_CUBE] = ExtractAverageCUBEPSFs(PSF_CUBE);

[PSF_poked_CUBE] = ExtractAverageCUBEPSFs(PSF_poked_CUBE);

% Average Together Taken Images, Center/Crop PSFs, Compute OTFs and dOTF
[dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF(PSF_CUBE, PSF_poked_CUBE,PT{1},PT{2});




end
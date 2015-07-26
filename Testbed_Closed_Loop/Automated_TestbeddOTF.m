function [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in)
% function [dOTF, PSF, PSF_poked, OTF, OTF_poked] = TestbeddOTF(DM, pokeact, verbose, Ppos_in)
%[dOTF, PSF, PSF_poked, OTF, OTF_poked] = TestbeddOTF(DM, pokeact, verbose, Ppos_in)
%
% Function to compute dOTF using the UAWFC testbed to obtain PSFs
% 
% INPUTS:
% DM is an AODM model of a 1024 pin BMC Mirror
% pokeact is the actuator to poke for the pupil blocker
% verbose is t/f flag to enable/disable plotting images
% Ppos_in is the initial actuator piston list
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
elseif nargin > 3
    error('Number of Inputs is Wrong: dOTF = TestbeddOTF(DM, pokeact,verbose, Ppos_in)');
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


%SEND TO MIRROR
% input('Press Enter to Send to Mirror');
fprintf('Sending Commands to DM\n');
! ~/src/scripts/dmloadch DM_flat.fits 0
cd(tempdir);


%CHECK FOR MIRROR UPDATE
% input('Press Enter once Mirror has Updated to Start Taking Images');
fprintf('Starting Camera...\n');
pause(2);
PSF_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);


% input('Press Enter to Modify Pupil');
fprintf('Modifying the Pupil\n');
% pause(1);
% ! ~/src/scripts/dmzeroch 0

%% Modified Setup
% Set DM and Take the Second PSF Image
Ppos_poked = Ppos_flat;
Ppos_poked(17,7) = Ppos_poked(17,7) + (AOField.HeNe_Laser*10^6) / 4;


cd /home/lab/src/scripts
fitswrite(Ppos_poked,'DM_poked.fits');


%SEND TO MIRROR
% input('Press Enter to Send to Mirror');
fprintf('Sending Commands to the DM\n');
! ~/src/scripts/dmloadch DM_poked.fits 0
cd(tempdir);

%CHECK FOR MIRROR UPDATE
% input('Press Enter once Mirror has Updated to Start Taking Images');
fprintf('Starting Camera...\n');
pause(2);
PSF_poked_CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters);

! ~/src/scripts/dmzeroch 0
%% Process the Pictures
[PSF_CUBE] = ExtractAverageCUBEPSFs(PSF_CUBE);

[PSF_poked_CUBE] = ExtractAverageCUBEPSFs(PSF_poked_CUBE);

% Average Together Taken Images, Center/Crop PSFs, Compute OTFs and dOTF
[dOTF, PSF_CUBE, PSF_poked_CUBE] = ExtractdOTF(PSF_CUBE, PSF_poked_CUBE);


% DM.setActs(Ppos_flat);


end
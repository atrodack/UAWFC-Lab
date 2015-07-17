function [dOTF, PSF_CUBE, PSF_poked_CUBE] = TestbeddOTF(DM,pokeact,verbose,Ppos_in)
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

if nargin < 4
    verbose = false;
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
elseif nargin < 5
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
else
    error('Number of Inputs is Wrong: dOTF = TestbeddOTF(DM, pokeact,verbose, Ppos_in)');
end

if ~isa(DM,'AODM')
    error('DM MUST BE an AODM Model of a 1024 Pin BMC Mirror: Class Type is Incorrect');
else
    if DM.nActs ~= 1024
        error('DM MUST BE an AODM Model of a 1024 Pin BMC Mirror: Number of Actuators in Model is Incorrect');
    end
end

if pokeact <= 1 || pokeact == 32 || pokeact == 993 || pokeact >= 1024
    error('pokeact MUST BE an Actuator on the Mirror: 2-1023, excluding 32 and 993');
end

if length(Ppos_in) ~= 1024
    error('Ppos_in MUST HAVE 1024 Elements');
end


%% Parameter Initialization
% Initializing a Couple Variables
varargin = cell(1,4);
Ppos_flat = Ppos_in;
Run_Cam_Parameters{1} = 10;
Run_Cam_Parameters{2} = 5;
Run_Cam_Parameters{3} = true;
Run_Cam_Parameters{4} = 1;
Run_Cam_Parameters{5} = 'HeNe';
Run_Cam_Parameters{6} = 'OD3';

%% Initial Setup
% Set DM and Take the First PSF Image
DM.setActs(Ppos_flat);
DM_Pistons = reshape(DM.actuators(:,3),[32,32]);
DM_Pistons(1,1) = 0; DM_Pistons(32,32) = 0; DM_Pistons(32,1) = 0; DM_Pistons(1,32) = 0;
DM_Pistons = single(DM_Pistons);


tempdir = pwd;
cd /home/alex/Desktop/Testbed_fits_files;
fitswrite(DM_Pistons,'DM_Pistons.fits');
img = fitsread('DM_Pistons.fits');

if verbose == true
    figure(5);
    imagesc(img);
    input('Press Enter');
    close
end

cd(tempdir);

% DO STUFF TO SEND TO MIRROR
input('Press Enter to Send to Mirror');


%CHECK FOR MIRROR UPDATE
input('Press Enter once Mirror has Updated to Start Taking Images');
PSF_CUBE = Testbed_run_cam(Run_Cam_Parameters);


input('Press Enter to Modify Pupil');

%% Modified Setup
% Set DM and Take the Second PSF Image
Ppos_poked = Ppos_flat;
Ppos_poked(pokeact) = Ppos_poked(pokeact) + (AOField.HeNe_Laser*10^6) / 4;

DM.setActs(Ppos_poked);
DM_Pistons_poked = reshape(DM.actuators(:,3),[32,32]);
DM_Pistons_poked(1,1) = 0; DM_Pistons_poked(32,32) = 0; DM_Pistons_poked(32,1) = 0; DM_Pistons_poked(1,32) = 0;
DM_Pistons_poked = single(DM_Pistons_poked);



cd /home/alex/Desktop/Testbed_fits_files;
fitswrite(DM_Pistons_poked,'DM_Pistons_poked.fits');
img = fitsread('DM_Pistons_poked.fits');

if verbose == true
    figure(5);
    imagesc(img);
    input('Press Enter');
    close
end

cd(tempdir);

% DO STUFF TO SEND TO MIRROR
input('Press Enter to Send to Mirror');


%CHECK FOR MIRROR UPDATE
input('Press Enter once Mirror has Updated to Start Taking Images');

PSF_poked_CUBE = Testbed_run_cam(Run_Cam_Parameters);


%% Process the Pictures
[PSF_CUBE] = ExtractAverageCUBEPSFs(PSF_CUBE);

[PSF_poked_CUBE] = ExtractAverageCUBEPSFs(PSF_poked_CUBE);

% Average Together Taken Images, Center/Crop PSFs, Compute OTFs and dOTF
dOTF = ExtractdOTF(PSF_CUBE, PSF_poked_CUBE);


DM.setActs(Ppos_flat);


end
% Closed-Loop dOTF Demo

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



%% Set Wavelength
% lambda = AOField.HeNe_Laser; 
lambda = 600 *10^-9;
bandpassfilter = 10; %0=no filter, 40=40nm band, 10=10nm band

% Compute Wavenumber
k = (2*pi) / lambda;

%% Set Camera Parameters
Run_Cam_Parameters{1} = 10; %Number of Triggers (0 based, want 10 triggers, this equals 9)
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

%% Data Loop
nn = 1;
niterations = 1;
grating_name = 1;
verbose = true;


while(nn <= niterations)
    fprintf('\nLoop Number %d\n',nn);
    
    % Measure the dOTF
    [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
    
    % Plot if asked to
    if verbose == true
        % Zero high amplitude pixels
        dOTF(129,:) = 0;
        dOTF = -1i * conj(dOTF);
        
        figure(10);
        plotComplex(dOTF,6);
        axis xy;
    end
    %% Write a Filename using the Date
    dt = datestr(now,'mm_dd_yyyy_HH_MM');
    filename = sprintf('Grating_%d_Testbed_Data_%d_%s',grating_name,nn,dt);
    % filename_movie_avi = sprintf('%s.avi',filename);
    fprintf('\n');
    display(filename);
    fprintf('\n');
    
    %% Save Cubes
    mkdir(filename);
    current_dir = pwd;
    cd(filename)
    savename = sprintf('Grating_%d_run_%d_Data_Cubes.mat',grating_name,nn);
    save(savename,'PSF_CUBE','PSF_poked_CUBE')
    cd(current_dir);

    
end



%% Plot Cubes
% figure(3);
% for n = 1:nn-1
%     imagesc(log10(normalize(psf_CUBE(:,:,n))),[-1.0 0]);
%     axis xy; sqar; colormap((gray));
%     drawnow;
%     pause(0.5);
% end
% 
% for n = 1:nn-1
%     imagesc(mirrorshape_CUBE(:,:,n));
%     sqar; colormap((gray));
%     drawnow;
%     pause(0.5);
% end
% % 
% for n = 1:nn-1
%     plotComplex(dOTF_CUBE(:,:,n),6);
%     axis xy; sqar; colormap((gray));
%     drawnow;
%     pause(0.5);
% end


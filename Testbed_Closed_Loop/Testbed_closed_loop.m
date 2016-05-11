clear all;
clc;
closepreview;
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



Run_Cam_Parameters{1} = 5;
Run_Cam_Parameters{2} = 30;
Run_Cam_Parameters{3} = 20;
Run_Cam_Parameters{4} = false;
Run_Cam_Parameters{5} = 1;
Run_Cam_Parameters{6} = 'SuperK';
Run_Cam_Parameters{7} = 'bp10nm';
Run_Cam_Parameters{8} = 16;
Run_Cam_Parameters{9} = 60;

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
% Ppos_in = fitsread('Testpattern1.fits');
% Ppos_in = fitsread('Testpattern1.fits');
% 
% 
% % Do a dOTF with the testpattern
% [dOTF_cal, PSF_CUBE_cal, PSF_poked_CUBE_cal] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in);
% 
% dOTF = dOTF_cal;
% dOTF(129,:) = 0;
% 
% phase = angle(dOTF_cal);
% uphase = uwrap(phase,'unwt');
% OPL = uphase / k;



% Plot
% figure;
% subplot(1,3,1);
% imagesc(PSF_CUBE_cal.PSF_centered_and_cropped); axis xy; axis off; sqar; bigtitle('PSF',15); %colormap(gray);
% 
% subplot(1,3,2)
% plotComplex(dOTF,10); axis xy; axis off; sqar; bigtitle('dOTF',15);
% 
% subplot(1,3,3);
% imagesc(OPL); axis xy; axis off; sqar; bigtitle('OPL',15);


% figure;
% plotComplex(dOTF,6);
% axis xy;
% sqar;
% % 
% % For Testpattern1, Pick center point, then point to right, then CCW
% PT = pickPoint(9);
% PT_map = [529 817 727 538 343 241 331 520 715];


%% dOTF Masks
[X,Y] = meshgrid(1:1:256);
R = sqrt((X-129).^2+ (Y-103).^2);
mask1 = double(R<=30);
R2 = sqrt((X-129).^2+ (Y-155).^2);
mask2 = double(R2<= 30);
mask = mask1 + mask2;
mask(mask~=0) = 1;





%% Closed Loop
nn = 1;
niterations = 5;
dOTF_CUBE = zeros(256,256,niterations);
psf_CUBE = zeros(256,256,niterations);
mirrorshape_CUBE = zeros(32,32,niterations);
% mirror_shape = (fitsread('Pupil_mask_Correction.fits'));
mirror_shape = zeros(32,32);
gain = -0.5;
load PT.mat;


%movie setup
% moviefig = figure(2);
% clf;
% input('Press Enter when figure is sized to liking');
% winsize = get(moviefig,'Position');
% winsize(1:2) = [0 0];
% MOVIE = moviein(niterations,moviefig,winsize);
% set(moviefig,'NextPlot','replacechildren');
% display('Movie Calibration Complete');
% fprintf('\n');

while(nn <= niterations)
    fprintf('\nLoop Number %d\n',nn);
    Ppos_in = mirror_shape;
    
    
    
    if nn == 1
        [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in,0);
        %     centerpoint = PSF_CUBE.centerpoint;
        %     croppoint = PSF_CUBE.croppoint;
        %     PT{1} = centerpoint;
        %     PT{2} = croppoint;
        
    else
        Run_Cam_Parameters{3} = 15;
        [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF(DM,Run_Cam_Parameters,Ppos_in,0);
    end
    
    dOTF(129,:) = 0;
    dOTF = -1i * conj(dOTF);
    
    figure(10);
    plotComplex(dOTF,6);
    axis xy;
    
    y = 'y';
    Y = 'y';
    n = 'n';
    N = 'n';
    display('Is this a usable dOTF?');
    user_check = input('y/n: ');
    if strcmpi(user_check,'y') == 1;
        dOTF_CUBE(:,:,nn) = dOTF;
        psf_CUBE(:,:,nn) = PSF_CUBE.PSF_centered_and_cropped;
        
        [ mirror_shape_new ] = computeAct_positions( dOTF );
        
        mirror_shape = (Ppos_in + (gain * mirror_shape_new));
        
        SELECT = mirror_shape(mirror_shape~=0);
        MEAN = mean(SELECT);
        for m = 1:32
            for n = 1:32
                if mirror_shape(m,n) ~= 0
                    mirror_shape(m,n) = mirror_shape(m,n) - MEAN;
                end
            end
        end
        mirrorshape_CUBE(:,:,nn) = mirror_shape;
        
        figure(2);
        subplot(1,3,1)
        imagesc(PSF_CUBE.PSF_centered_and_cropped);
        axis xy;
        sqar;
        title('PSF');
        
        subplot(1,3,2)
        plotComplex(dOTF,6);
        % imagesc(OPL *1e6);
        axis xy;
        sqar;
        title('dOTF');
        
        subplot(1,3,3)
        imagesc(mirror_shape);
        axis ij; axis off;
        sqar;
        title('Shape to be Sent to Mirror');
        drawnow;
        
%         figure(2)
%         subplot(1,3,1)
%         imagesc(abs(dOTF));
%         axis xy;
%         sqar;
%         
%         subplot(1,3,2);
%         imagesc(angle(dOTF));
%         axis xy;
%         sqar;
%         
%         subplot(1,3,3);
%         % imagesc(mirror_shape_new);
%         plotComplex(dOTF,6);
%         axis xy;
%         sqar;
%         % colorbar;
%         % caxis([-0.5 0.5])
        
        fprintf('Maximum Magnitude in mirror_shape is: %0.5f \n',max(max(abs(mirror_shape))));
        % pause(30);
        % input('Continue?');
        % MOVIE(:,nn) = getframe(moviefig,winsize);
        nn = nn+1
        
    elseif strcmpi(user_check,'n') == 1
        fprintf('Starting loop over to recompute dOTF\n\n');
    end
    
end




mkdir(filename);
current_dir = pwd;
cd(filename)
save('Data_Cubes.mat','PSF_CUBE','PSF_poked_CUBE','dOTF_CUBE','psf_CUBE','mirrorshape_CUBE');
cd(current_dir);

% %%
% figure(3);
% for n = 1:nn-1
%     imagesc(log10(normalize(psf_CUBE(:,:,n))),[-1.5 0]);
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


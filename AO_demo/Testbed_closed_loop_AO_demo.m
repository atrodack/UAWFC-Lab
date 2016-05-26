% Closed-Loop dOTF Demo

%% System Parameters

% Set Wavelength
% lambda = AOField.HeNe_Laser; 
lambda = 600 *10^-9;
bandpassfilter = 10; %0=no filter, 40=40nm band, 10=10nm band

% Compute Wavenumber
k = (2*pi) / lambda;

% Set Camera Parameters
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

%% Closed Loop
nn = 1;
niterations = 10;
dOTF_CUBE = zeros(256,256,niterations);
psf_CUBE = zeros(256,256,niterations);
mirrorshape_CUBE = zeros(32,32,niterations);
mirror_shape = zeros(32,32);
gain = -0.75;

while(nn <= niterations)
    fprintf('\nLoop Number %d\n',nn);
    Ppos_in = mirror_shape;
    
    % Measure the dOTF (and lower the shutter speed after first correction)
    if nn == 1
        [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
    elseif nn < 6
        Run_Cam_Parameters{3} = 20;
        [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
    else
        Run_Cam_Parameters{3} = 10;
        [dOTF, PSF_CUBE, PSF_poked_CUBE] = Automated_TestbeddOTF_AO_demo(DM,Run_Cam_Parameters,Ppos_in,0,[17,8],(lambda*10^6)/2);
    end
    
    % Zero high amplitude pixels
    dOTF(129,:) = 0;
    dOTF = -1i * conj(dOTF);
    
    figure(10);
    plotComplex(dOTF,6);
    axis xy;
    
    y = 'y';
    Y = 'y';
    n = 'n';
    N = 'n';
%     display('Is this a usable dOTF?');
%     user_check = input('y/n: ');
    user_check = 'y';
    if strcmpi(user_check,'y') == 1;
        dOTF_CUBE(:,:,nn) = dOTF;
        psf_CUBE(:,:,nn) = PSF_CUBE.PSF_centered_and_cropped;
        
        [ mirror_shape_new ] = computeAct_positions_AO_demo( dOTF, k );
        
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
        clf;
        title('dOTF was used')
        axis off
        
        figure(6);
        subplot(1,3,1)
        imagesc(PSF_CUBE.PSF_centered_and_cropped);
        axis xy;
        sqar;
        title('PSF');
        
        subplot(1,3,2)
        plotComplex(dOTF,2);
        % imagesc(OPL *1e6);
        axis xy;
        sqar;
        title(sprintf('dOTF for iteration %d',nn));
        
        subplot(1,3,3)
        imagesc(mirror_shape);
        axis ij; axis off;
        sqar;
        title('Shape to be Sent to Mirror');
        drawnow;
        
        
        fprintf('Maximum Magnitude in mirror_shape is: %0.5f \n',max(max(abs(mirror_shape))));
        nn = nn+1;
        
    elseif strcmpi(user_check,'n') == 1
        fprintf('Starting loop over to recompute dOTF\n\n');
    end
    
end

%% Save Cubes
mkdir(filename);
current_dir = pwd;
cd(filename)
save('Data_Cubes.mat','PSF_CUBE','PSF_poked_CUBE','dOTF_CUBE','psf_CUBE','mirrorshape_CUBE');
cd(current_dir);

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


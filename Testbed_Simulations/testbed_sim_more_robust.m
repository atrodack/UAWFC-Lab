clear all;
clc;
close all;

%**************************************************************************
%                       UAWFC Testbed Simulation
%**************************************************************************

%% System Parameters
global lambda k D secondary spider SPACING aa fftsize THld FOV PLATE_SCALE FoV_withIrisAO FoV_withoutIrisAO RunSIM RunTESTBED IrisAO_on BMC_on verbose_makeDM Scalloped_Field UseRealPSF coronagraph system_verbose
% lambda = AOField.RBAND; % Red light.
lambda = AOField.HeNe_Laser;
k = (2*pi) / lambda;
D = 7e-3; % 7mm
secondary = 0.3*D; % 30% of Pupil Diameter
spider = 0.02*D; % 2% of Pupil Diameter

%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 4; %number of zernikes to inject
goal_strehl = 0.9; %exit condition for loop
gain = 0.7;
fftsize = 2^12;

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 100*THld; % arcsecs
PLATE_SCALE = THld/4;
FoV_withIrisAO = 3.5e-3;
FoV_withoutIrisAO = 10.5e-3;

%% Run Simulation/Testbed Flags
RunSIM = true; %Run the simulation
RunTESTBED = false; %Run the testbed equipment

%% Simulation Flags

% IrisAO Flags
IrisAO_on = false;
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.
% BMC Flag
BMC_on = true;
% Aberration Flag
InjectAb = true; %Injects some random Zernikes
InjectRandAb = true;
InjectKnownAb = false;
% Corrector Flag
UseDM4Correction = false;
% Noise Flags
UseNoise = false;
if UseNoise == true
    number_of_images = 100;
    Noise = cell(2,1);
    Noise{1} = UseNoise;
    Noise{2} = number_of_images;
else
    Noise = cell(1,1);
    Noise{1} = UseNoise;
end
% Use Testbed PSF Instead of Simulated PSF
UseRealPSF = true;
if UseRealPSF == true
    InjectAb = false;
    Num_Folders = 2;
    Num_files_per_folder = 100;
%     varargin{1} = '/home/alex/Desktop/Data/2015612_Batch1_nofilter_PSFWithoutFinger/';
%     varargin{3} = 'RAW_scienceIM_frame_';
%     varargin{2} = '/home/alex/Desktop/Data/2015612_Batch2_nofilter_PSFWithFinger/';
%     varargin{4} = 'RAW_scienceIM_frame_';
    varargin{1} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithoutFingerDMBox/';
    varargin{3} = 'RAW_scienceIM_frame_';
    varargin{2} = '/home/alex/Desktop/Data/2015615_Batch1_nofilter_PSFWithFingerDMBox/';
    varargin{4} = 'RAW_scienceIM_frame_';
end
% Coronagraph Flag
coronagraph = false; % turns on going through coronagraph elemens
% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Flags



%% Make the Testbed Elements
maketestbedelements;



%% Initialize dOTF Object
if RunSIM == true 
    if RunTESTBED == true
        dOTF_Sim = AOdOTF(A,FOV,PLATE_SCALE,FoV,SPACING);
        dOTF_Sim.lambda = lambda;
        dOTF_Testbed = AOdOTF(1);
        dOTF_Testbed.lambda = lambda;
        fprintf('\ndOTF Objects Created for Sim and Testbed\n');
        dOTF_Sim.create_finger(pupil_blocker(1),pupil_blocker(2),pupil_blocker(3));
    else
        dOTF_Sim = AOdOTF(A,FOV,PLATE_SCALE,FoV,SPACING);
        dOTF_Sim.lambda = lambda;
        fprintf('\ndOTF Object Created for Sim\n');
        dOTF_Sim.create_finger(pupil_blocker(1),pupil_blocker(2),pupil_blocker(3));
        
    end
else
    if RunTESTBED == true
        dOTF_Testbed = AOdOTF(1);
        dOTF_Testbed.lambda = lambda;
        fprintf('\ndOTF Object Created for Testbed\n');
    else
        fprintf('WHY WOULD YOU EVEN BOTHER RUNNING THIS?\n');
    end
end

if UseRealPSF == true
    dOTF_Sim.usedata = true;
end

%% Inject Zernike Aberration to simulate system errors
if RunSIM == true
    if InjectAb == true && InjectRandAb == true
        ABER = AOScreen(A);
        n = sort((randi(4,1,nzerns)),'ascend'); %zernike mode order (from lowest to highest)
        m = zeros(1,nzerns); %zernike mode initialization
        
        %Get a correct "m" index for each "n" index
        for ii = 1:nzerns
            m_pos = -n(ii):2:n(ii);
            choice = randi(length(m_pos),1,1);
            m(ii) = m_pos(choice);
        end
        
        % Find a coefficient between -1 and 1 waves, scale it to the number of
        % zernikes added to avoid unrealistic PSF deaths (fairly arbitrary,
        % probably uncessary)
        if nzerns < 3
            coeffs = (2*rand(1,nzerns)-1);
        elseif nzerns <= 5
            coeffs = (((nzerns-(nzerns-2))/nzerns)) .* (2*rand(1,nzerns)-1);
        elseif nzerns > 5
            coeffs = (((nzerns-(nzerns-5))/nzerns)) .* (2*rand(1,nzerns)-1);
        end
        
        % Add the Zernikes into ABER
        ABER.zero;
        for ii = 1:nzerns
            ABER.addZernike(n(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n = n';
        m = m';
        Number_of_waves = coeffs';
        T = table(n,m,Number_of_waves);
        fprintf('\nInjected Aberrations:\n');
        disp(T);
    elseif InjectAb == true && InjectKnownAb == true
        ABER = AOScreen(A);
        n = [2,2,2,3,3];
        m = [-2,0,2,-1,3];
        
%         coeffs = 1 * randn(1,length(n));
        coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
        ABER.zero;
        for ii = 1:length(n)
            ABER.addZernike(n(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n = n';
        m = m';
        Number_of_waves = coeffs';
        T = table(n,m,Number_of_waves);
        fprintf('\nInjected Aberrations:\n');
        disp(T);
    end
end

%% Get the Diffraction Limited PSF and Calibrate WFS
F = AOField(A);
F.spacing(SPACING);
F.FFTSize = fftsize; % Used to compute PSFs, etc.
F.resize(F.FFTSize);
F.planewave;
F.lambda = lambda;
F * A * DM1 * DM2;

% dOTF_Sim.calibrateWFS2(F);
dOTF_Sim.precalibratedWFS(1);


DLPSF = F.mkPSF(FOV,PLATE_SCALE);

DLmax = max(DLPSF(:));
clear F;

%% Make the Coronagraph Elements
if RunSIM == true
    if coronagraph == true
        F_coronagraph = AOField(A);
        F_coronagraph.FFTSize = fftsize;
        F_coronagraph.spacing(SPACING);
        F_coronagraph.planewave * A;
        F_coronagraph.grid(F_coronagraph.fft/F_coronagraph.nx); % Go to the focal plane.
        
        FPMASK = AOSegment(F_coronagraph);
        FPMASK.spacing(SPACING);
%         FPMASK.grid(exp(-normalize(F_coronagraph.mag2)/0.05) ); % This is pretty ad hoc.
        
        %Do it "Better"
        D_FPM = 10e-5;
        
        FPM_DEFN = [
        0 0 D_FPM         1 aa 0 0 0 0 0
        ];
    
        FPMASK.pupils= FPM_DEFN;
        FPMASK.make;
        fpmask = FPMASK.grid;
        fpmask = padarray(fpmask,[floor((fftsize-length(fpmask))/2),floor((fftsize-length(fpmask))/2)]);
        fpmask = double(~fpmask);
        FPMASK.grid(fpmask);
%         fprintf('\nFPM Constructed\n\n');
        
        F_coronagraph.grid(F_coronagraph.fft/F_coronagraph.nx); % Go to the Lyot pupil plane.
        
        LYOT = AOSegment(F_coronagraph);
        LYOTSTOP_DEFN = [
            0 0 (D*0.8)         1 aa 0 0 0 0 0  % undersize the Lyot stop
            0 0 (secondary*1.1) 0 aa/2 0 0 0 0 0 % oversize the secondary
            0 0 spider         -2 aa 4 0 D/1.9 0 0
            ];
        
        LYOT.pupils = LYOTSTOP_DEFN;
        LYOT.make;
%         fprintf('Lyot Stop Constructed\n');
        
        [PSF0,thx,thy] = F_coronagraph.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
        PSFmax = max(PSF0(:)); % Save for normalizing.
        
        PSF0 = PSF0/PSFmax; % make the brightest value =1.
        fprintf('Coronagraph Elements Constructed\n');
    end
end
%% *************************************************************************
%                    Model Light through the System
%**************************************************************************
fprintf('\nSending Light through the System and computing dOTF\n');
n = 1;
strehl = 0;
figure(1);
drawnow;
CORRECTOR = 1;
while(n<=5)
    fprintf('\nLoop # %d\n',n);
    if RunSIM == true
        if UseRealPSF == false
            if Scalloped_Field == true
                fprintf('Using a Scalloped Field to Simulate Segment Surfaces\n');
                F = F_scal.copy;
                F.FFTSize = fftsize; % Used to compute PSFs, etc.
            else
                F = AOField(A);
                F.spacing(SPACING);
                F.FFTSize = fftsize; % Used to compute PSFs, etc.
                F.resize(F.FFTSize);
                F.planewave;
            end
            F.lambda = lambda;
            
            
            if InjectAb == true
                if UseDM4Correction == true
                    F * ABER * A * DM1 * DM2;
                else
                    F * ABER * A * DM1 * CORRECTOR;
                end
            else
                if UseDM4Correction == true
                    F * A * DM1 * DM2;
                else
                    F * A * DM1 * CORRECTOR;
                end
            end
            
            figure(1);
            subplot(2,3,1);
            F.show;
            title('Field Through System');
            
            if coronagraph == true
                dOTF_Sim.sense_coronagraph(F,FPMASK,LYOT); %doesn't seem to work right for no IrisAO, but dOTF shouldn't be done with coronagraph in anyway
            else
                dOTF_Sim.sense2(F,Noise,'unwt');
            end
            
            if n == 1
                PSF = dOTF_Sim.PSF0;
                PSFmax2 = max(PSF(:));
            end
            
            
            if n > 10
                strehl_previous = strehl;
            end
            
            strehl = (max(max(dOTF_Sim.PSF0))/DLmax);
            fprintf('Tip/Tilt Insensitive Strehl is %0.5f\n',strehl);
        else
            dOTF_Sim.useData2(Num_Folders,Num_files_per_folder,false,varargin);
            PSF = dOTF_Sim.PSF0;
            PSFmax2 = max(PSF(:));
            mag = abs(dOTF_Sim.dOTF);
            mag(129,129) = 0;
            
            subplot(1,3,1);
            imagesc(PSF);
            axis off;
            sqar;
            title('Testbed PSF');
            
            subplot(1,3,2);
            imagesc(mag);
            axis off;
            sqar;
            title('dOTF Magnitude');
            
            subplot(1,3,3)
            imagesc(dOTF_Sim.Phase);
            axis off;
            sqar;
            title('dOTF Phase');
            
        end
        
        
        %% Test the dOTF Result
        if UseRealPSF == false
            CORRECTOR = AOScreen(1);
            CORRECTOR.spacing(SPACING);
            OPL = dOTF_Sim.OPL;
            OPL(OPL~=0) = OPL(OPL~=0)-mean(mean(OPL));
            OPL = padarray(OPL,[floor((length(DM2.grid)-length(OPL))/2),floor((length(DM2.grid)-length(OPL))/2)],'both');
            CORRECTOR.grid(OPL);
%             CORRECTOR * A;
            pistonvec = CORRECTOR.interpGrid(DM2.actuators(DM2.OnActs,1),DM2.actuators(DM2.OnActs,2));
            DM2.bumpOnActs(gain*pistonvec);
            storeDMcommands{n} = DM2.actuators(:,3);
            DM2.clip(STROKE);
            DM2.removeMean;
            DM2.render;
            
            
            
            subplot(2,3,2);
            surf(x_DM,y_DM,DM2.grid,'LineStyle','none');
            zlim([-1.5e-6,1.5e-6]);
            daspect([1 1 1e-3]);
            lt=light();
            set(lt,'Position',[-15e-3 0 15e-3]);
            colormap(gray);
            
            subplot(2,3,3)
            imagesc(x,y,dOTF_Sim.OPL);
            axis xy;
            colormap(gray);
            %         caxis([-1.5e-6,1.5e-6]);
            sqar;
            title('dOTF Computed OPL');
            
            subplot(2,2,3)
            imagesc(dOTF_Sim.thx,dOTF_Sim.thy,log10(PSF/DLmax),[-4,0]);
            %         imagesc(dOTF_Sim.thx,dOTF_Sim.thy,log10(PSF));
            axis xy;
            colormap(gray);
            sqar;
            title('Uncorrected PSF');
            
            subplot(2,2,4)
            imagesc(dOTF_Sim.thx,dOTF_Sim.thy,log10(dOTF_Sim.PSF0 / DLmax),[-4,0]);
            %         imagesc(dOTF_Sim.thx,dOTF_Sim.thy,log10(dOTF_Sim.PSF0));
            axis xy;
            colormap(gray);
            sqar;
            title(sprintf('Loop # %d PSF',n));
            drawnow;
            
                        
            n = n+1;
            
        else
            n = 1000;
            fprintf('****** Make sure to pick Bottom Pupil Center First *******\n');
            dOTF_Sim.mkMask;
            phase_copy = dOTF_Sim.Phase;
            phase = dOTF_Sim.Phase;
            mask = dOTF_Sim.Mask;
            
            shift_point = dOTF_Sim.pupil_center;
            radius = round(dOTF_Sim.pupil_radius);
            [sizex,sizey] = size(phase);
            center = [round(sizex/2),round(sizey/2)];
            shift = [shift_point(1) - center(2),shift_point(2) - center(1)];
            phase_ref = phase(center(1),center(2));
            
            if dOTF_Sim.calibration == false
                fprintf('\nUnwrapping the Phase\n');
                dOTF_Sim.unwrapphase('unwt');
                phase = dOTF_Sim.Phase - phase_ref;
            end
            
            phase = circshift(phase,-shift);
            phase_copy = circshift(phase_copy,-shift);
            
            if dOTF_Sim.calibration == false
                fprintf('Cropping, Masking, and Resizing Computed Phase\n');
            else
                fprintf('Creating Correctly Scaled Mask\n');
            end
            
            % resize to edge of Pupil
            phase = phase(center(1)-radius-0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);
            phase_copy = phase_copy(center(1)-radius-0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);


            dOTF_Sim.Phase = phase;
            dOTF_Sim.resize_phase_to_Pupil;
            phase = dOTF_Sim.Phase .* dOTF_Sim.Mask_interped;
            dOTF_Sim.calibration = true;
            dOTF_Sim.Phase = phase_copy;
            dOTF_Sim.resize_phase_to_Pupil;
            phase_copy = dOTF_Sim.Phase .* dOTF_Sim.Mask_interped;
            dOTF_Sim.Phase = phase;
            
            
            fprintf('Computing the OPL\n');
            dOTF_Sim.OPL = dOTF_Sim.Phase / k;
            
            [Ax,Ay] = A.coords;
            figure;
            imagesc(Ax,Ay,dOTF_Sim.OPL);
            sqar;
            colorbar;
            axis xy;
%             DM2.plotActuators(1);
            
            
            
            CORRECTOR = AOScreen(1);
            CORRECTOR.spacing(SPACING);
            OPL = dOTF_Sim.OPL;
            OPL(OPL~=0) = OPL(OPL~=0)-mean(mean(OPL));
            OPL = padarray(OPL,[floor((length(DM2.grid)-length(OPL))/2),floor((length(DM2.grid)-length(OPL))/2)],'both');
            CORRECTOR.grid(OPL);
%             CORRECTOR * A;
            pistonvec = CORRECTOR.interpGrid(DM2.actuators(DM2.OnActs,1),DM2.actuators(DM2.OnActs,2));
            DM2.bumpOnActs(gain*pistonvec);
            storeDMcommands{n} = DM2.actuators(:,3);
            DM2.clip(STROKE);
            DM2.removeMean;
            DM2.render;
            figure;
            DM2.show;
%             DM2.plotActuators(1);
            drawnow;
%             dOTF_Sim.scanNumericalDefocus(-1.5,0.5,1000);
        end
    end
end



load ok.mat;
John = audioplayer(y,Fs);
% play(John);
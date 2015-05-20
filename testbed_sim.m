clear all;
clc;
close all;

%**************************************************************************
%                       UAWFC Testbed Simulation
%**************************************************************************

%% System Parameters
lambda = AOField.RBAND; % Red light.
D = 7e-3; % 7mm
secondary = 0.3*D; % 30% of Pupil Diameter
spider = 0.02*D; % 2% of Pupil Diameter

%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 9; %number of zernikes to inject
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
% Aberration Flag
InjectAb = true; %Injects some random Zernikes
% Coronagraph Flag
coronagraph = false; % turns on going through coronagraph elemens
% Plotting Flag
system_verbose = false; %Plots Created System Elements
% Testbed Flags


%% ************************************************************************
%                      Construct System Elements
%**************************************************************************
if RunSIM == true
    load alright.mat;
    John = audioplayer(y,Fs);
    play(John);
    
    %% Pupil Mask
    PUPIL_DEFN = [
        0 0 D         1 aa 0 0 0 0 0
        0 0 secondary 0 aa/2 0 0 0 0 0
        0 0 spider   -2 aa 4 0 D/1.9 0 0
        ];
    
    A = AOSegment;
    A.spacing(SPACING);
    A.name = 'UAWFC Pupil Mask';
    A.pupils = PUPIL_DEFN;
    A.make;
    [x,y] = A.coords;
    fprintf('\nPupil Mask Constructed\n');
    
    if system_verbose == true
        figure(1)
        A.show;
        colormap(gray);
        drawnow;
        pause(1);
    end
end

%% IrisAO DM
if RunSIM == true
    if IrisAO_on == true
        segpitch = 606e-6; %leave this alone
        magnification = 1; %leave this alone too
        FoV = FoV_withIrisAO;
        % Make the Mirror
        if Scalloped_Field == true
            [DM1,F_scal] = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
            F_scal.grid(padarray(F_scal.grid,[ceil(271/2),ceil(207/2)]));
        else
            DM1 = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
        end
        DM1.lambdaRef = lambda;
        
        % clc; %clears the command window of the text generated from adding segments to the AOAperture object
        fprintf('\nIrisAO Mirror Constructed\n');
        pupil_blocker = [1.788e-3,0.3212e-3,0.1e-3];
    else
        DM1 = 1;
        FoV = FoV_withoutIrisAO;
        Scalloped_Field = false;
        fprintf('\nUsing a Flat instead of the IrisAO\n');
        pupil_blocker = [3.28e-3,0.3212e-3,0.1e-3];
    end
end

%% Set the Initial Piston, Tip, Tilt of IrisAO Mirror
% Flatten the IrisAO
PTTpos = zeros(37,3);

% Set a Random Mirror Shape
% PTTpos = horzcat(randn(37,1)*10^-6,randn(37,1)*10^-3,randn(37,1)*10^-3);

% Use a Zernike on the mirror....magnification must be equal to 1 for this!
% Zernike_Number = [2,3];
% Zernike_Coefficient_waves = randn(2,1);
% PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);

if RunTESTBED == true
    fprintf('*******************************\nSending PTTpos to IrisAO Mirror\n*******************************\n\n');
    %     tempdir = pwd;
    %     cd /home/lab/Desktop/Shared_Stuff
    %     save('PTTpos','PTTpos');
    %     checkmirrorupdated = false;
    %     while(checkmirrorupdated == false)
    %         CMD_FILES = dir('PTTpos.mat');
    %         if(~isempty(CMD_FILES))
    %             pause(0.1);
    %         else
    %             checkmirrorupdated = true;
    %         end
    %     end
    %     cd(tempdir);
    %     clear tempdir;
end

if RunSIM == true
    if IrisAO_on == true
        % Load in Mapping Data
        load('IrisAO_SegMap.mat');
        PTT = zeros(37,3);
        
        % Map the PTT matrix from hardware to software order
        for ii = 1:37
            mapped_segment = IrisAO_SegMap(ii);
            PTT(ii,1:3) = PTTpos(mapped_segment,:);
        end
        
        % Send to DM Model
        DM1.PTT(PTT);
        DM1.touch;
        DM1.render;
        
        if system_verbose == true
            figure(1);
            DM1.show;
            colormap(gray);
            drawnow;
            pause(1);
        end
    end
end

%% Boston MircroMachines DM
%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
STROKE = 1.5e-6;
Pitch = 300e-6;
radius_BMC = (4.65)*10^-3;

if RunSIM == true
    %Coordinate System
    xmin = -radius_BMC;
    xmax = radius_BMC;
    BMC_x = (xmin:SPACING:xmax);
    ymin = -radius_BMC;
    ymax = radius_BMC;
    BMC_y = (ymin:SPACING:ymax);
    [BMC_X,BMC_Y] = meshgrid(BMC_x,BMC_y);
    
    BMC_pupil = ones(size(BMC_X));
    % BMC_pupil = padarray(BMC_pupil,[250,250],1,'both'); %put 250 pixels on both sides outside of active area
    
    %Construct the Pupil
    Seg = AOSegment(length(BMC_pupil));
    Seg.spacing(SPACING);
    Seg.name = 'BMC Pupil';
    Seg.grid(BMC_pupil);
    
    %Make it an AOAperture Class
    A_BMC = AOAperture;
    A_BMC.spacing(SPACING);
    A_BMC.name = 'BMC Aperture';
    A_BMC.addSegment(Seg);
    
    %Make an AODM
    DM2 = AODM(A_BMC);
    [X,Y] = DM2.COORDS;
    
    % Create Actuator Locations
    actuator_x = xmin:Pitch:xmax;
    actuator_x = actuator_x - mean(actuator_x);
    actuator_y = ymin:Pitch:ymax;
    actuator_y = actuator_y - mean(actuator_y);
    [ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
%     ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
%     ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
%     ACTUATOR_X(ACTUATOR_X==0) = [];
%     ACTUATOR_Y(ACTUATOR_Y==0) = [];
    BMC_ACTS = [ACTUATOR_X(:),ACTUATOR_Y(:)];
    nacts = length(BMC_ACTS(:,1));
    
    % Add Actuators
    DM2.addActs(BMC_ACTS,1,A_BMC);
    
    % Turn Off Actuators that Aren't Illuminated
    RHO = zeros(nacts,1);
    for ii = 1:nacts
        RHO(ii) = sqrt(BMC_ACTS(ii,1)^2 + BMC_ACTS(ii,2)^2);
        if RHO(ii) > D/1.8
            DM2.actuators(ii,5) = 0;
        elseif RHO(ii) < secondary/2.1
%             DM2.actuators(ii,5) = 0;
        end
    end
    
    % Get List of Which Actuators are Being Used
    DM2.setOnActs;
    
    %Turn Off Actuators so the Program Knows they are off
    DM2.disableActuators(DM2.OffActs);
    
    % Set the Convex Hull Boundary Conditions
    DM2.defineBC(D/2.1,108,'circle');
    
    [x_DM,y_DM] = DM2.coords;
    %% Set the Initial Piston Values of the BMC DM
    DM2.flatten;
%     DM2.actuators(343,3) = 1e-6;
%     DM2.actuators(OnActs,3) = ((randn(length(OnActs),1).*10^-6));
    DM2.removeMean;
    
    if RunTESTBED == true
       DM_Pistons = reshape(DM2.actuators(:,3),[32,32]);
       
       tempdir = pwd;
       cd C:\Users\atrod_000\Documents\GitHub\UAWFC\fits_files;
       fitswrite(DM_Pistons,'DM_Pistons.fits');
       img = fitsread('DM_Pistons.fits');
       figure(5);
       imagesc(img);
       input('Press Enter');
       close
       cd(tempdir);
       
       % DO STUFF TO SEND TO MIRROR
       
       
    end
    
    if system_verbose == true
        DM2.show;
        colormap(gray);
        title('BMC Pupil');
        DM2.plotActuators;
        pause(1);
        DM2.plotRegions;
        pause(1);
    end
    
    fprintf('\nBMC DM Constructed\n');
end

%% Initialize dOTF Object
if RunSIM == true 
    if RunTESTBED == true
        dOTF_Sim = AOdOTF(A,FOV,PLATE_SCALE,FoV,SPACING);
        dOTF_Testbed = AOdOTF(1);
        fprintf('\ndOTF Objects Created for Sim and Testbed\n');
        dOTF_Sim.create_finger(pupil_blocker(1),pupil_blocker(2),pupil_blocker(3));
    else
        dOTF_Sim = AOdOTF(A,FOV,PLATE_SCALE,FoV,SPACING);
        fprintf('\ndOTF Object Created for Sim\n');
        dOTF_Sim.create_finger(pupil_blocker(1),pupil_blocker(2),pupil_blocker(3));
        
    end
else
    if RunTESTBED == true
        dOTF_Testbed = AOdOTF(1);
        fprintf('\ndOTF Object Created for Testbed\n');
    else
        fprintf('WHY WOULD YOU EVEN BOTHER RUNNING THIS?\n');
    end
end

%% Inject Random Zernike Aberration to simulate system errors
if RunSIM == true
    if InjectAb == true
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
while(strehl < goal_strehl)
    fprintf('\nLoop # %d\n',n);
    if RunSIM == true
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
            F * ABER * A * DM1 * DM2;
        else
            F * A * DM1 * DM2;
        end
        
        figure(1);
        subplot(2,3,1);
        F.show;
        title('Field Through System');
        
        if coronagraph == true
            dOTF_Sim.sense_coronagraph(F,FPMASK,LYOT); %doesn't seem to work right for no IrisAO, but dOTF shouldn't be done with coronagraph in anyway
        else
            dOTF_Sim.sense2(F,'gold');
        end
        
        if n == 1
            PSF = dOTF_Sim.PSF0;
            PSFmax2 = max(PSF(:));
        end
        
%         if system_verbose == true
%             PSF1 = dOTF_Sim.PSF1;
%             PSFmax3 = max(PSF1(:));
%             thx = dOTF_Sim.thx;
%             thy = dOTF_Sim.thy;
%             OTF = dOTF_Sim.OTF0;
%             OTF1 = dOTF_Sim.OTF1;
%             OTFmax = max(abs(OTF(:)));
%             OTF1max = max(abs(OTF1(:)));
%             dOTF = dOTF_Sim.dOTF;
%             
%             figure(2)
%             subplot(1,2,1)
%             imagesc(thx,thy,PSF/PSFmax2);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('PSF\n'));
%             
%             subplot(1,2,2)
%             imagesc(thx,thy,log10(PSF/PSFmax2),[-3,0]);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('Log Scale PSF\n'));
%             input('Press Enter');
%             
%             figure(2)
%             clf;
%             subplot(1,2,1)
%             plotCAmpl(OTF/OTFmax,1);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('OTF\n'));
%             
%             subplot(1,2,2)
%             plotCAmpl(log10(OTF/OTFmax),0.5);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('Log Scale OTF\n'));
%             input('Press Enter');
%             
%             figure(3)
%             subplot(1,2,1)
%             imagesc(thx,thy,PSF1/PSFmax3);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('PSF\n'));
%             subplot(1,2,2)
%             imagesc(thx,thy,log10(PSF1/PSFmax3),[-3,0]);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('Log Scale PSF\n'));
%             input('Press Enter');
%             
%             figure(3)
%             clf;
%             subplot(1,2,1)
%             plotCAmpl((OTF1/OTF1max),1);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('OTF\n'));
%             
%             subplot(1,2,2)
%             plotCAmpl(log10(OTF1/OTF1max),0.5);
%             colormap(gray);
%             axis xy;
%             sqar;
%             title(sprintf('Log Scale OTF\n'));
%         end
        
        if n > 10
            strehl_previous = strehl;
        end
            
        strehl = (max(max(dOTF_Sim.PSF0))/DLmax);
        fprintf('Tip/Tilt Insensitive Strehl is %0.5f\n',strehl);
        
        
%         if strehl > goal_strehl
%             fprintf('The Strehl is above %0.2f, breaking loop\n',goal_strehl);
%             strehl_final = strehl;
%             break;
%         end
        
        %% Test the dOTF Result
        CORRECTOR = AOScreen(1);
        CORRECTOR.spacing(SPACING);
        OPL = dOTF_Sim.OPL;
        OPL(OPL~=0) = OPL(OPL~=0)-mean(mean(OPL));
        OPL = padarray(OPL,[floor((length(DM2.grid)-length(OPL))/2),floor((length(DM2.grid)-length(OPL))/2)],'both');
        CORRECTOR.grid(OPL);
%         CORRECTOR * A;
        pistonvec = CORRECTOR.interpGrid(DM2.actuators(DM2.OnActs,1),DM2.actuators(DM2.OnActs,2));
        DM2.bumpOnActs(gain*pistonvec);
        storeDMcommands{n} = DM2.actuators(:,3);
        DM2.clip(1.5*STROKE);
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
        axis xy;
        colormap(gray);
        sqar;
        title('Uncorrected PSF');
        
        subplot(2,2,4)
        imagesc(dOTF_Sim.thx,dOTF_Sim.thy,log10(dOTF_Sim.PSF0 / DLmax),[-4,0]);
        axis xy;
        colormap(gray);
        sqar;
        title(sprintf('Loop # %d PSF',n));
        drawnow;
        
        if n > 20
            if(strehl - strehl_previous < 0.001)
                strehl_final = strehl;
                strehl = 1;
            end
        end
        
        n = n+1;
    end
end
if n>20
    strehl = strehl_final;
end


load ok.mat;
John = audioplayer(y,Fs);
play(John);
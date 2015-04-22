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
aa = SPACING;  % for antialiasing.

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 300; % arcsecs
PLATE_SCALE = THld/5;

%% Run Simulation/Testbed Flags
RunSIM = true; %Run the simulation
RunTESTBED = false;%Run the testbed equipment

%% Simulation Flags
% IrisAO Flags
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.                  
% Aberration Flag
InjectAb = true; %Injects some random Zernikes
% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Flags


%**************************************************************************
%                      Construct System Elements
%**************************************************************************
if RunSIM == true
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
    segpitch = 606e-6; %leave this alone
    magnification = 1; %leave this alone too

    % Make the Mirror
    if Scalloped_Field == true
        [IrisAO_DM,F_scal] = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
    else
        IrisAO_DM = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
    end
    IrisAO_DM.lambdaRef = lambda;
    
    % clc; %clears the command window of the text generated from adding segments to the AOAperture object
    fprintf('\nIrisAO Mirror Constructed\n');
end
%% Set the Initial Piston, Tip, Tilt of IrisAO Mirror

% Flatten the IrisAO
PTTpos = zeros(37,3);

% Set a Random Mirror Shape
% PTTpos = horzcat(randn(37,1)*10^-6,randn(37,1)*10^-3,randn(37,1)*10^-3);

% Use a Zernike on the mirror
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
    % Load in Mapping Data
    load('IrisAO_SegMap.mat');
    PTT = zeros(37,3);
    
    % Map the PTT matrix from hardware to software order
    for ii = 1:37
        mapped_segment = IrisAO_SegMap(ii);
        PTT(ii,1:3) = PTTpos(mapped_segment,:);
    end
    
    % Send to DM Model
    IrisAO_DM.PTT(PTT);
    IrisAO_DM.touch;
    IrisAO_DM.render;
    
    if system_verbose == true
        figure(1);
        IrisAO_DM.show;
        colormap(gray);
        drawnow;
        pause(1);
    end
end

%% Boston MircroMachines DM
%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
Max_Stroke = 1.5e-6;
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
    BMC_DM = AODM(A_BMC);
    [X,Y] = BMC_DM.COORDS;
    
    % Create Actuator Locations
    actuator_x = xmin:Pitch:xmax;
    actuator_y = ymin:Pitch:ymax;
    [ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
    ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
    ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
    ACTUATOR_X(ACTUATOR_X==0) = [];
    ACTUATOR_Y(ACTUATOR_Y==0) = [];
    BMC_ACTS = zeros(1020,2);
    BMC_ACTS(:,1) = ACTUATOR_X(:);
    BMC_ACTS(:,2) = ACTUATOR_Y(:);
    
    % Add Actuators
    BMC_DM.addActs(BMC_ACTS);
    BMC_DM.defineBC(radius_BMC+Pitch,7,'square');
    BMC_DM.flatten;
    
    if system_verbose == true
        BMC_DM.show;
        colormap(gray);
        title('BMC Pupil');
        BMC_DM.plotActuators;
        pause(1);
        BMC_DM.plotRegions;
        pause(1);
    end
    
    fprintf('\nBMC DM Constructed\n');
end

%% Initialize dOTF Object
if RunSIM == true 
    if RunTESTBED == true
        dOTF_Sim = AOdOTF(A_BMC,FOV,PLATE_SCALE);
        dOTF_Testbed = AOdOTF(1);
        fprintf('\ndOTF Objects Created for Sim and Testbed\n');
    else
        dOTF_Sim = AOdOTF(A_BMC,FOV,PLATE_SCALE);
        fprintf('\ndOTF Object Created for Sim\n');
    end
else
    if RunTESTBED == true
        dOTF_Testbed = AOdOTF(1);
        fprintf('\ndOTF Object Created for Testbed\n');
    else
        while(1)
            fprintf('WHY WOULD YOU EVEN BOTHER RUNNING THIS?\n');
            pause(0.15);
        end
    end
end


%% Inject Random Zernike Aberration to simulate system errors
if RunSIM == true
    if InjectAb == true
        nzerns = 2; %number of zernikes to add
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
        if nzerns == 1
            coeffs = (2*rand(1,nzerns)-1);
        elseif nzerns > 1
            coeffs = (2/nzerns) .* (2*rand(1,nzerns)-1);
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
%**************************************************************************
%                    Model Light through the System
%**************************************************************************
if RunSIM == true
    if Scalloped_Field == true
        F = F_scal.copy;
        F.FFTSize = 2048; % Used to compute PSFs, etc.
    else
        F = AOField(A);
        F.spacing(SPACING)
        F.FFTSize = 2048; % Used to compute PSFs, etc.
        F.resize(F.FFTSize);
        F.planewave;
    end
    F.lambda = lambda;
    
    if InjectAb == true
        F * ABER * A * IrisAO_DM * BMC_DM;
    else
        F * A * IrisAO_DM * BMC_DM;
    end
    
    figure(1);
    F.show;
    drawnow;
    title('Field After IrisAO DM');
    
    figure(2)
    [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
    PSFmax = max(PSF(:));
    subplot(1,2,1)
    imagesc(thx,thy,PSF/PSFmax);
    colormap(gray);
    axis xy;
    sqar;
    title(sprintf('PSF\n'));
    subplot(1,2,2)
    imagesc(thx,thy,log10(PSF/PSFmax),[-2,0]);
    colormap(gray);
    axis xy;
    sqar;
    title(sprintf('Log Scale PSF\n'));
    
end
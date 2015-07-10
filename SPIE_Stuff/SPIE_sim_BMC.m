clearvars -except ActMasks
clc;
% close all;

% Closed-Loop dOTF Methods and Testing

%**************************************************************************
%                       SPIE dOTF Closed-Loop Methods Simulation
%**************************************************************************

%% System Parameters
% Don't touch these if using maketestbedelements.m to construct components
global lambda k D secondary spider SPACING aa fftsize THld FOV PLATE_SCALE FoV_withIrisAO FoV_withoutIrisAO RunSIM RunTESTBED IrisAO_on BMC_on verbose_makeDM Scalloped_Field UseRealPSF coronagraph system_verbose

% Set Wavelength
% lambda = AOField.RBAND; % Red light.
lambda = AOField.HeNe_Laser;

% Compute Wavenumber
k = (2*pi) / lambda;

% Pupil Specs
D = 7e-3; % 7mm
% secondary = 0.3*D; % 30% of Pupil Diameter
% spider = 0.02*D; % 2% of Pupil Diameter



%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 9; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
numiterations = 100;
gain = 0.7; %gain for AO Corrections
fftsize = 2^11;

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 35*THld; % arcsecs
PLATE_SCALE = THld/4;
FoV_withIrisAO = 3.5e-3;
FoV_withoutIrisAO = 10.5e-3;

%% Run Simulation/Testbed Flags
RunSIM = true; %Run the simulation
RunTESTBED = false; %Run the testbed equipment

%% Simulation Flags

% IrisAO Flags
IrisAO_on = false; %turns on/off IrisAO Mirror (if false, DM1 variable set to 1)
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.

% BMC Flag
BMC_on = true; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = false; %Injects nzerns Zernike Terms
InjectRandAb = false; %if InjectAB is true, picks Zernikes "Randomly"
InjectKnownAb = true; %if InjectAB is true, picks provided Zernikes

InjectKolm = true;

% Check Aberration Flags
if InjectKolm == true
    if InjectAb == true
        y = 'y';
        Y = 'y';
        n = 'n';
        N = 'n';
        display('You have turned on both Zernikes and Kolmogorov Turbulence. Would you like to turn one off?');
        answer1 = input('y/n : ');
        if strcmpi(answer1,'y') == true
            Zern = 'Zern';
            zern = 'Zern';
            Kolm = 'Kolm';
            kolm = 'Kolm';
            display('Turn off which Aberration Type?');
            answer2 = input('Zern/Kolm : ');
            if strcmpi(answer2,'Zern') == true
                InjectAb = false;
                display('Zernikes will no longer be injected');
            elseif strcmpi(answer2,'Kolm') == true
                InjectKolm = false;
                display('Kolmogorov Turbulence will no longer be injected');
            else
                error('Incorrect Input: Must be either Zern or Kolm');
            end
            clear Zern Kolm zern kolm
        end
        clear y n Y N
    end
end

% Corrector Flag
UseDM4Correction = true;

% Noise Flags
UseNoise = false;
if UseNoise == true
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 5;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 0.5;
    Noise_Parameters{4} = 0.5;
    Noise_Parameters{5} = UseNoise;
else
    Noise_Parameters = cell(5,1);
%     Noise_Parameters{1} = 1;
%     Noise_Parameters{2} = false;
%     Noise_Parameters{3} = 0;
%     Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
end

% Use Testbed PSF Instead of Simulated PSF
UseRealPSF = false;
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

% Plotting Flag
system_verbose = false; %Plots Created System Elements

%% Testbed Flags



%% Make the Testbed Elements
maketestbedelements;


%% Inject Aberration
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
        % zernikes added to avoid unrealistic PSFs (fairly arbitrary,
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
        %         n = [2,2,2,3,3];
        %         n = [1,1,2,4];
        n = [1];
        %         m = [-2,0,2,-1,3];
        %         m = [-1,1,0,0];
        m = [-1];
        
        %         coeffs = 1 * randn(1,length(n));
        %         coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
        %         coeffs = 0.25*randn(1,length(n));
        coeffs = [4];
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
        %         figure(1);
        %         ABER.show;
        
        
    elseif InjectAb == false
        ABER = 1;
    end
    
    if InjectKolm == true
        TURB = AOScreen(A,4e-3,lambda);
        TURB.spacing(SPACING);
        TURB.name = 'Simulated Turbulence';
        TURB.make;
        wind_dir = randn(2,1);
        wind_dir = wind_dir./abs(wind_dir);
        wind_strength = randi(10,2,1);
        Wind = wind_dir .* wind_strength
        
        r0 = TURB.r0;
        D_fit = 7e-3;
        %  expected_strehl = exp(-(1.03 * (D_fit/r0) ^ (5/3))) %not tip/tilt corrected
        expected_strehl = exp(-(0.13 * (D_fit/r0) ^ (5/3)))
    else
        TURB = 1;
        
    end
end

%% BMC Simulation

% load in calibration file
% load BMC_calibrated_X_and_Y_dOTF_locations
% calibrated_BMC_act_locations = cell(DM2.nActs,1);
% counter = 1;
% for n = 1:32
%     for m = 1:32
%         hold on
%         calibrated_BMC_act_locations{counter} = [X(m,n),Y(m,n)];
%         drawnow;
%         pause(0.1);
%         hold off
%         counter = counter + 1;
%     end
% end
load calibrated_BMC_dOTF_act_locations;
for n = 1:DM2.nActs
    if DM2.actuators(n,5) == 0
        calibrated_BMC_act_locations{n} = [];
    end
end

onAct_locations =  calibrated_BMC_act_locations(~cellfun('isempty',calibrated_BMC_act_locations));


% make a field object
display('Making Fields');
new_spacing = DM2.spacing;
F = AOField(fftsize);
F.FFTSize = (fftsize);
F.spacing(new_spacing);
F.lambda = lambda;
F.FOV = FOV;
F.PLATE_SCALE = PLATE_SCALE;
F.FoV = FoV_withoutIrisAO;
F.name = 'IrisAO Field 1';

% make diffraction limited PSF
display('Making Diffraction Limited PSF');
F.planewave * A * DM2;
[PSF_difflim,plotx,ploty] = F.mkPSF(FOV,PLATE_SCALE);
PSF_difflimmax = max(max(PSF_difflim));
F.touch;

strehl_location = ceil(size(PSF_difflim)/2);


% DM2.flatten;
% DM2.poke(497,3e-6);
% ABER.grid(DM2.grid);


%% Control Loop
nn = 1;
no_actuators = [1,32,993,1024];
% DM2.actuators(no_actuators,5) = 0;
Ppos = zeros(DM2.nActs,1);

% sizey = fftsize;
% sizex = fftsize;
% xx = linspace(1,sizex,sizex);
% yy = linspace(1,sizey,sizey);
% [XX,YY] = meshgrid(xx,yy);
% 
% ActMasks = cell(length(DM2.OnActs),1);
% for n = 1:length(DM2.OnActs)
%         R = sqrt((XX - onAct_locations{n}(2)).^2 + (YY - onAct_locations{n}(1)).^2);
%         ActMasks{n,1} = double(R <= 5);
% end

figure(1);
drawnow;
input('Press Enter when figure is sized to liking');

Ppos2 = zeros(length(onAct_locations),1);
% load ActMasks
while(nn <= numiterations)
    
    if InjectKolm == true
        TURB.grid(circshift(TURB.grid,Wind));
    end
    Ppos = DM2.actuators(:,3);
    DM2.flatten;
    F.planewave * TURB *ABER * A * DM2;
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;
    
    [ dOTF, PSF1, PSF2, OTF1, OTF2 ] = BMCcomputedOTF( DM2, 698, Ppos, Noise_Parameters, F, A, ABER, TURB );
    
    mag = abs(dOTF);
    maxmag = max(max(mag));
    thresh = maxmag / 1000;
    mask = mag;
    mask(mask<thresh) = 0;
    mask = double(mask>0);
    
    G = AOGrid(1);
    G.spacing(DM2.spacing);
    G.grid(dOTF);
    [X, Y] = G.COORDS;
    GT = G.copy;
    GT.grid(~(Y - tand(150)*X > -0.0004)); % fftsize = 2^10, manual ffts
%     GT.show
    G * GT;
    
    
    
    masked_dOTF = G.grid;
    dOTF = -1i * conj(dOTF);
    
    phase = angle(dOTF);
    unwrapped_phase = uwrap(phase,'gold');
    OPL = unwrapped_phase / k;
    OPL = OPL .* GT.grid;

    

    
%     for n = 1:length(DM2.OnActs)
%         piston_area = OPL .* ActMasks{n};
%         Ppos2(n,1) = mean(mean(abs(piston_area)>0));
%     end
    
    for n = 1:length(DM2.OnActs)
        Ppos2(n,1) = OPL(onAct_locations{n}(2),onAct_locations{n}(1));
    end
    
    SELECT = Ppos2(Ppos2~=0);
    MEAN = mean(SELECT);
    for n = 1:length(Ppos2)
        if Ppos2(n) ~= 0
            Ppos2(n) = Ppos2(n) - MEAN;
        end
    end
    
    %% Correction
    if nn > 1  %suffer the seeing limit for a bit
        gain = -0.7;
        DM2.bumpOnActs(gain * (Ppos2));
        
        DM2.render;
    else
        
    end
    
    F.touch;
    F.planewave * TURB *ABER * A * DM2;
    % figure(2);
    % F.show;
    
    PSF_cor = F.mkPSF(FOV,PLATE_SCALE);
    maxPSF_cor = max(max(PSF_cor));
    
    
    [DM2x,DM2y] = DM2.coords;
    
    figure(1);
    subplot(2,4,1)
    imagesc(plotx,ploty,log10(PSF_difflim / PSF_difflimmax),[-4,0]);
    sqar;
    title('Diffraction Limited PSF');
    subplot(2,4,2)
    imagesc(plotx,ploty,log10(PSF_aberrated / PSF_aberratedmax),[-4,0])
    sqar;
    title('Aberrated PSF')
    subplot(2,4,5);
    imagesc(plotx,ploty,log10(PSF_cor / maxPSF_cor),[-4,0]);
    sqar;
    title(sprintf('Corrected PSF, loop #%d',nn));
    subplot(2,4,6);
    strehl(nn) = PSF_cor(strehl_location(1),strehl_location(2)) / PSF_difflim(strehl_location(1),strehl_location(2));
    strehl_uncorr(nn) = PSF_aberrated(strehl_location(1),strehl_location(2)) / PSF_difflim(strehl_location(1),strehl_location(2));
    strehl_notip(nn) = maxPSF_cor / PSF_difflimmax;
    loopnum(nn) = nn;
    plot(loopnum,strehl,'-r');
    hold on
    plot(loopnum,strehl_uncorr,'-b');
    plot(loopnum,strehl_notip,'-g');
    hold off
    xlabel('Loop Iteration');
    ylabel('Strehl Ratio');
    legend('Strehl for Corrected PSF','Strehl for Uncorrected/Aberrated PSF','Tip/Tilt Independent Strehl','Location','Best');
    xlim([0,100]);
    ylim([0,1]);
    title('Strehl Ratio');

    subplot(2,4,3);
    A_mask = A.copy;
    A_mask * DM2;
    A_mask.show;
    colorbar;
    axis xy;
    sqar;
    bigtitle(sprintf('DM Shape with Correction Applied\n'),10);
    subplot(2,4,4);
    F.planewave * TURB * ABER * A * DM2;
    F.show;
    axis xy;
    sqar;
    colorbar off;
    bigtitle('Field at DM',10);
    subplot(2,4,7)
    plotComplex(masked_dOTF,5);
    axis off;
    axis xy;
    sqar;
    colorbar off;
    bigtitle(sprintf('dOTF and Segment Center Locations\n'),10);
%     hold on
%     for n = 1:length(DM2.OnActs)
%         plot(onAct_locations{n}(1),onAct_locations{n}(2),'g.');
%     end
%     hold off
    subplot(2,4,8)
    if InjectKolm == true
        TURB.show;
        sqar;
        colorbar;
        bigtitle(sprintf('Turbulence Profile at loop %d \n',nn),10);
        xlim([-2e-3,2e-3]);
        ylim([-2e-3,2e-3]);
    else
        ABER.show;
        colorbar;
        bigtitle(sprintf('Injected Aberration at loop %d \n',nn),10);
        xlim([-4e-3,4e-3]);
        ylim([-4e-3,4e-3]);
        sqar;
    end
%     
%     
    drawnow;
    
%     figure(2);
%     subplot(1,2,1)
%     DM2.show;
%     subplot(1,2,2)
%     imagesc(unwrapped_phase);
%     drawnow;
% %     
%     fprintf('Approximate Strehl: %0.5f \n',strehl(nn));
    nn = nn + 1;
end


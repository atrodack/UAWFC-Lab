clear all;
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
D = 5e-3; % 7mm
% secondary  0.3*D; % 30% of Pupil Diameter
secondary = 0;
% spider = 0.02*D; % 2% of Pupil Diameter
spider = 0;


%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 4; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
% numiterations = 2;
gain = 0.7; %gain for AO Corrections
fftsize = 2^11;

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
IrisAO_on = true; %turns on/off IrisAO Mirror (if false, DM1 variable set to 1)
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.

% BMC Flag
BMC_on = false; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = false; %Injects nzerns Zernike Terms
InjectRandAb = true; %if InjectAB is true, picks Zernikes "Randomly"
InjectKnownAb = false; %if InjectAB is true, picks provided Zernikes

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
UseNoise = true;
if UseNoise == true
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 5;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 0.75;
    Noise_Parameters{4} = 0;
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
% Coronagraph Flag
coronagraph = false; % turns on going through coronagraph elemens

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
%                 n = [2,2,2,3,3];
        %         n = [1,1,2,4];
        n = [1];
%                 m = [-2,0,2,-1,3];
        %         m = [-1,1,0,0];
        m = [1];
        
%                 coeffs = 1 * randn(1,length(n));
        %         coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
%                 coeffs = 0.25*randn(1,length(n));
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
        wobble_dir = randn(2,1);
        wobble_dir = wobble_dir./abs(wobble_dir);
        wobble_strength = randi(3,2,1);
        Wobble = wobble_dir .* wobble_strength;
        
    elseif InjectAb == false
        ABER = 1;
    end
    
    if InjectKolm == true
        
        %         TURB = AOAtmo(A);
        % %         TURB.spacing(SPACING);
        %         WFlow = AOScreen(fftsize,0.15,500e-9);
        %         WFlow.spacing(SPACING);
        %         WFlow.name = 'Lower altitude turbulence';
        %         WFhigh = AOScreen(2*fftsize,0.17,500e-9);
        %         WFhigh.spacing(SPACING);
        %         WFhigh.name = 'High altitude turbulence';
        %
        %         TURB.addLayer(WFlow,1000);
        %         TURB.addLayer(WFhigh,8000);
        %
        %         TURB.layers{1}.Wind = [3 1];
        %         TURB.layers{2}.Wind = [1 -1]*20;
        %
        %         r0 = TURB.totalFriedScale;
        %         th_scat = lambda/r0*206265;
        %
        %         fprintf('The total r0 is %f cm.\n',100*r0);
        %         fprintf('The seeing is %.2f arcsecs.\n',th_scat);
        %
        %         % Turning this off is like using dynamic refocus.
        %         TURB.GEOMETRY = false;
        %         TURB.BEACON = [1 1 1e10];
        %         TURB.make;
        % %         TURB.show
        %
        %         aberration = TURB.grid;
        
        TURB = AOScreen(A,4e-3,lambda);
        TURB.spacing(SPACING);
        TURB.name = 'Simulated Turbulence';
        TURB.make;
%         grid = TURB.grid;
%         TURB.grid(grid * 30);
        wind_dir = randn(2,1);
        wind_dir = wind_dir./abs(wind_dir);
        wind_strength = randi(5,2,1);
        Wind = wind_dir .* wind_strength
        
        r0 = TURB.r0;
        D_fit = 4e-3;
        %  expected_strehl = exp(-(1.03 * (D_fit/r0) ^ (5/3)))
        expected_strehl = exp(-(0.13 * (D_fit/r0) ^ (5/3)))
    else
        TURB = 1;
        
    end
end

%% IrisAO Simulation

% [pixel_seg_map,Areal_Averaging_radius] = computeIrisAOsegpixelmap(DM1, A, 23, FOV, PLATE_SCALE, FoV_withIrisAO);
calibration_filename_23 = 'calibrated_dOTF_segment_centers_for_phase_finger_segment_23_no_aliasing.mat';
calibration_filename_35 = 'calibrated_dOTF_segment_centers_for_phase_finger_segment_35_no_aliasing.mat';
load(calibration_filename_23);
% load(calibration_filename_35);

new_spacing = DM1.spacing;

display('Making Fields');

F = AOField(fftsize);
F.FFTSize = (fftsize);
F.spacing(new_spacing);
F.lambda = lambda;
F.FOV = FOV;
F.PLATE_SCALE = PLATE_SCALE;
F.FoV = FoV_withIrisAO;
F.name = 'IrisAO Field 1';

display('Making Diffraction Limited PSF');
F.planewave * A * DM1;
[PSF_difflim,plotx,ploty] = F.mkPSF(FOV,PLATE_SCALE);
PSF_difflimmax = max(max(PSF_difflim));
F.touch;


PTTpos_mirror = zeros(37,3);
PTT_flat = zeros(37,3);
%% Control Loop
nn = 1;
DM1.isMirror = 0;
numiterations = 100;
%movie setup
moviefig = figure(1);
input('Press Enter when figure is sized to liking');
winsize = get(moviefig,'Position');
winsize(1:2) = [0 0];
MOVIE = moviein(numiterations,moviefig,winsize);
set(moviefig,'NextPlot','replacechildren');






while(nn <= numiterations)
    %     if InjectAb == true
    %         wobble_dir = randn(2,1);
    %         wobble_dir = wobble_dir./abs(wobble_dir);
    %         wobble_strength = randi(3,2,1);
    %         Wobble = wobble_dir .* wobble_strength;
    %         ABER.grid(circshift(ABER.grid,Wobble));
    %     end
    
    if InjectKolm == true
        TURB.grid(circshift(TURB.grid,Wind));
    end
    
    DM1.setIrisAO(PTT_flat);
    
    F.planewave * ABER * TURB * A * DM1;
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;
    
    [ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 23, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
    dOTF = -1i * conj(dOTF);
    
    mag = abs(dOTF);
    %     phase = angle(dOTF);
    %     unwrapped_phase = uwrap(phase,'unwt');
    %     OPL = unwrapped_phase / k;
    
    maxmag = max(max(mag));
    thresh = maxmag / 10;
    mask = mag;
    mask(mask<thresh) = 0;
    mask = double(mask>0);
    
    
    
    
    
    %% Correction
    if nn > 1  %suffer the seeing limit for a bit
        
        [ PTT_mirror, PTTpos_mirror1] = IrisAOgetPTT_new( dOTF, [1,1], lambda, [22,23,24], 23, calibration_filename_23);
        
        %% Add Second dOTF to get overlap Region
        
        [ dOTF2, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 35, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
        dOTF2 = conj(dOTF2);
        [ PTT_mirror2, PTTpos_mirror2] = IrisAOgetPTT_new( dOTF2, [1,1], lambda, [34,35,36], 35, calibration_filename_35);
        %
        PTTpos_mirror3 = PTTpos_mirror1;
        PTTpos_mirror3([22,23,24],:) = PTTpos_mirror2([22,23,24],:);
        PTTpos_mirror3 = -PTTpos_mirror3;
        PTTpos_mirror3(:,1) = PTTpos_mirror3(:,1) * 0.654008264745021162;
        PTTpos_mirror = PTTpos_mirror + PTTpos_mirror3;
        PTT_mirror = mapSegments(PTTpos_mirror);
        DM1.setIrisAO(PTT_mirror);
    else
        PTTpos_mirror = zeros(37,3);
        DM1.setIrisAO(PTTpos_mirror);
    end
    
    F.touch;
    F.planewave * ABER * TURB * A * DM1;
    % figure(2);
    % F.show;
    
    PSF_cor = F.mkPSF(FOV,PLATE_SCALE);
    maxPSF_cor = max(max(PSF_cor));
    
    
    [DM1x,DM1y] = DM1.coords;
    
%     figure(1);
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
    strehl(nn) = PSF_cor(401,401) / PSF_difflim(401,401);
    strehl_uncorr(nn) = PSF_aberrated(401,401) / PSF_difflim(401,401);
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
    xlim([0,numiterations]);
    ylim([0,1]);
    title('Strehl Ratio');
    %     drawnow;
    %
    %
    %     figure(2)
    subplot(2,4,3);
    imagesc(DM1x,DM1y,angle(DM1.grid));
    colorbar;
    axis xy;
    sqar;
    bigtitle(sprintf('Phase of DM with Correction Applied\n'),10);
    subplot(2,4,4);
    F.planewave * TURB * ABER * DM1;
    F.show;
    axis xy;
    sqar;
    colorbar;
    bigtitle('Field at DM',10);
    subplot(2,4,7)
    plotComplex(dOTF,3);
    axis off;
    axis xy;
    sqar;
    colorbar;
    bigtitle(sprintf('dOTF and Segment Center Locations\n'),10);
    hold on
    for n = 1:37
        if n ~= 23
            plot(pixel_seg_map{n}(2),pixel_seg_map{n}(1),'r*');
        else
        end
    end
    hold off
    subplot(2,4,8)
    if InjectKolm == true
        TURB.show;
        sqar;
        bigtitle(sprintf('Turbulence Profile at loop %d \n',nn),10);
        xlim([-2e-3,2e-3]);
        ylim([-2e-3,2e-3]);
    else
        ABER.show;
        bigtitle(sprintf('Injected Aberration at loop %d \n',nn),10);
        xlim([-2e-3,2e-3]);
        ylim([-2e-3,2e-3]);
        sqar;
    end
    
    
    drawnow;
    MOVIE(:,nn) = getframe(moviefig,winsize);
    fprintf('Loop #%d Approximate Strehl: %0.6f \n',nn,strehl(nn));
    nn = nn + 1;
end


save closed_loop_movie.mat MOVIE winsize
% movie(figure(1),MOVIE,1,10,winsize)

% Zernike_Number = [2];
% Zernike_Coefficient_waves = 4;
% PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);
% PTT = mapSegments(PTTpos);
% pistonlist_zern = PTT(:,1);
% pistonlist = PTT_mirror(:,1);
% pistonlist = horzcat(pistonlist,pistonlist_zern);
% tiplist = horzcat(PTT_mirror(:,2),PTT(:,2));
% tiltlist = horzcat(PTT_mirror(:,3),PTT(:,3));

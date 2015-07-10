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
FOV = 50*THld; % arcsecs
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
InjectAb = true; %Injects nzerns Zernike Terms
InjectRandAb = false; %if InjectAB is true, picks Zernikes "Randomly"
InjectKnownAb = true; %if InjectAB is true, picks provided Zernikes

InjectKolm = false;

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

% Noise Flags
UseNoise = true;
if UseNoise == true
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 5;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 0.5;
    Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = 1.5e6;
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
        n_zern = sort((randi(4,1,nzerns)),'ascend'); %zernike mode order (from lowest to highest)
        m = zeros(1,nzerns); %zernike mode initialization
        
        %Get a correct "m" index for each "n" index
        for ii = 1:nzerns
            m_pos = -n_zern(ii):2:n_zern(ii);
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
            ABER.addZernike(n_zern(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n_zern = n_zern';
        m = m';
        Number_of_waves = coeffs';
        T = table(n_zern,m,Number_of_waves);
        fprintf('\nInjected Aberrations:\n');
        disp(T);
    elseif InjectAb == true && InjectKnownAb == true
        ABER = AOScreen(A);
        %                 n_zern = [2,2,2,3,3];
                n_zern = [1,1,2,4];
%         n_zern = [1];
        %                 m = [-2,0,2,-1,3];
                m = [-1,1,0,0];
%         m = [1];
        
        %  coeffs = 1 * randn(1,length(n_zern));
        %  coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
        %  coeffs = 0.25*randn(1,length(n_zern));
        coeffs = [0.0661 	0.1545 	-0.0022 	0.1995 	];
        %         coeffs = [4];
        ABER.zero;
        for ii = 1:length(n_zern)
            ABER.addZernike(n_zern(ii),m(ii),coeffs(ii)*lambda,D);
        end
        n_zern = n_zern';
        m = m';
        Number_of_waves = coeffs';
        T = table(n_zern,m,Number_of_waves);
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
    
    if InjectAb == true && InjectKolm == true
        PHASESCREEN = TURB.copy;
        TURB = 1;
    end
    
    
end



%% Write a Filename
if Noise_Parameters{5} == true
    s1 = sprintf('%d_photons',Noise_Parameters{6});
    s3 = sprintf('%0.2f_ReadNoise',Noise_Parameters{3});
    s4 = sprintf('%d_images_per_PSF',Noise_Parameters{1});
else
    s1 = 'No_Noise_Used';
    s3 = '_';
    s4 = '_';
end
if InjectAb == true
    s2 = sprintf('Zernike');
    if InjectKolm == true
        s2 = sprintf('Zernike_with_Kolmogorov_added');
    end
else
    if InjectKolm == true
        s2 = sprintf('Kolmogorov');
    else
        s2 = sprintf('No_Aberration');
    end
end


dt = datestr(now,'mm_dd_yyyy_HH_MM_SS_FFF');
filename = sprintf('Closed_loop_%s_%s_and_%s_with_%s',s2,s1,s3,s4);
filename = sprintf('%s_%s',filename,dt) %add date and time to msec in order to avoid overwriting files
filename_movie_save = sprintf('%s.mat',filename);
filename_movie_avi = sprintf('%s.avi',filename);
filename_txt = sprintf('%s.txt',filename);
% filename_OTFa = sprintf('%s_OTFaCUBE.mat',filename);
% filename_OTFb = sprintf('%s_OTFbCUBE.mat',filename);
filename_PTT = sprintf('%s_PTTCUBE.mat',filename);



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
fprintf('\n');
F.planewave * A * DM1;
[PSF_difflim,plotx,ploty] = F.mkPSF(FOV,PLATE_SCALE);
PSF_difflimmax = max(max(PSF_difflim));
F.touch;


PTTpos_mirror = zeros(37,3);
PTT_flat = zeros(37,3);
%% Control Loop
nn = 1;
DM1.isMirror = 0;
numiterations = 11;%60;
correction_start = 4;%6;
pixelshift = [1,1];

%movie setup
moviefig = figure(1);
clf;
input('Press Enter when figure is sized to liking');
winsize = get(moviefig,'Position');
winsize(1:2) = [0 0];
MOVIE = moviein(numiterations,moviefig,winsize);
set(moviefig,'NextPlot','replacechildren');
display('Movie Calibration Complete');
fprintf('\n');


%initalize data cubes and Strehl vectors
fprintf('Initializing Data Cubes\n\n');
% dOCUBE = zeros(fftsize,fftsize);
% PSFaCUBE = zeros(fftsize,fftsize);
% PSFbCUBE = 0;
% OTFaCUBE = zeros(fftsize,fftsize);
% OTFbCUBE = zeros(fftsize,fftsize);
DMCOMMANDSCUBE = zeros(37,3);

strehl = zeros(1,1);
strehl_notip = zeros(1,1);
strehl_uncorr = zeros(1,1);

load hexmask.mat; %removes smoothed edges from variance computation


fprintf('Starting the loop\n\n');
fprintf('Pixel Shift for slopes calculation: [%d,%d]\n\n',pixelshift(1),pixelshift(2));


while(nn <= numiterations)
    %     if InjectAb == true
    %         wobble_dir = randn(2,1);
    %         wobble_dir = wobble_dir./abs(wobble_dir);
    %         wobble_strength = randi(3,2,1);
    %         Wobble = wobble_dir .* wobble_strength;
    %         ABER.grid(circshift(ABER.grid,Wobble));
    %     end
    
    if InjectKolm == true
        if InjectAb == true
            PHASESCREEN.grid(circshift(PHASESCREEN.grid,Wind));
        else
            TURB.grid(circshift(TURB.grid,Wind));
        end
    end
    
    %Comput the uncorrected PSF
    DM1.setIrisAO(PTT_flat);
    F.planewave * ABER * TURB * A * DM1;
    W = angle(F.grid).*hexmask;
    WFE_uncor(nn) = var(W(abs(W)>0));
    strehl_uncorr(nn) = exp(-WFE_uncor(nn));
    
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;
    
    %Compute the dOTF
    [ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 23, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
    dOTF = -1i * conj(dOTF);
    
    
    %Store the dOTF and PSF1 data for access later
%     PSFaCUBE(:,:,nn) = PSF1;
%     PSFbCUBE(:,:,nn) = PSF2;
%     OTFaCUBE(:,:,nn) = OTF1;
%     OTFbCUBE(:,:,nn) = OTF2;
%     dOCUBE(:,:,nn) = dOTF;
    
    
    
    %% Correction
    if nn < correction_start  %suffer the seeing limit for a bit
        %Flatten the DM for the seeing limit
        PTTpos_mirror = zeros(37,3);
        PTT_mirror = zeros(37,3);
        DM1.setIrisAO(PTT_mirror);
        
    elseif nn >= correction_start && nn < 3*correction_start %close the loop
        %Get Segment PTT from dOTF
        [ ~, PTTpos_mirror1] = IrisAOgetPTT_new( dOTF, pixelshift, lambda, [22,23,24], 23, calibration_filename_23);
        
        %% Add Second dOTF to get overlap Region
        [ dOTF2, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 35, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
        dOTF2 = conj(dOTF2);
        [ PTT_mirror2, PTTpos_mirror2] = IrisAOgetPTT_new( dOTF2, pixelshift, lambda, [34,35,36], 35, calibration_filename_35);
        
        %Replace Overlap Region with data from second dOTF computation
        PTTpos_mirror3 = PTTpos_mirror1;
        PTTpos_mirror3([22,23,24],:) = PTTpos_mirror2([22,23,24],:);
        
        %Change the sign to correct
        PTTpos_mirror3 = -PTTpos_mirror3;
        
        %Apply the necessary scaling factor for the Piston
        PTTpos_mirror3(:,1) = PTTpos_mirror3(:,1) * 0.654008264745021162;
        
        %Bump the Segment PTT's
        PTTpos_mirror = PTTpos_mirror + PTTpos_mirror3;
        
        %Map them to code order
        PTT_mirror = mapSegments(PTTpos_mirror);
        
        %Apply the Correction to mirror
        DM1.setIrisAO(PTT_mirror);
        
    elseif nn >= 3*correction_start && nn < 6*correction_start
        if InjectAb == true && InjectKolm == true %Simulate adding in some hair spray
            TURB = PHASESCREEN.copy;
            [ ~, PTTpos_mirror1] = IrisAOgetPTT_new( dOTF, pixelshift, lambda, [22,23,24], 23, calibration_filename_23);
            [ dOTF2, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 35, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
            dOTF2 = conj(dOTF2);
            [ PTT_mirror2, PTTpos_mirror2] = IrisAOgetPTT_new( dOTF2, pixelshift, lambda, [34,35,36], 35, calibration_filename_35);
            PTTpos_mirror3 = PTTpos_mirror1;
            PTTpos_mirror3([22,23,24],:) = PTTpos_mirror2([22,23,24],:);
            PTTpos_mirror3 = -PTTpos_mirror3;
            PTTpos_mirror3(:,1) = PTTpos_mirror3(:,1) * 0.654008264745021162;
            PTTpos_mirror = PTTpos_mirror + PTTpos_mirror3;
            PTT_mirror = mapSegments(PTTpos_mirror);
            DM1.setIrisAO(PTT_mirror);
        else %Turn off correction for a bit
            %Flatten the DM
            PTTpos_mirror = zeros(37,3);
            PTT_mirror = zeros(37,3);
            DM1.setIrisAO(PTT_mirror);
        end

%         
    elseif nn >= 6*correction_start
        if InjectAb == true && InjectKolm == true %Simulate taking the hair spray out
            TURB = 1;
            
            [ ~, PTTpos_mirror1] = IrisAOgetPTT_new( dOTF, pixelshift, lambda, [22,23,24], 23, calibration_filename_23);
            [ dOTF2, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 35, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
            dOTF2 = conj(dOTF2);
            [ PTT_mirror2, PTTpos_mirror2] = IrisAOgetPTT_new( dOTF2, pixelshift, lambda, [34,35,36], 35, calibration_filename_35);
            PTTpos_mirror3 = PTTpos_mirror1;
            PTTpos_mirror3([22,23,24],:) = PTTpos_mirror2([22,23,24],:);
            PTTpos_mirror3 = -PTTpos_mirror3;
            PTTpos_mirror3(:,1) = PTTpos_mirror3(:,1) * 0.654008264745021162;
            PTTpos_mirror = PTTpos_mirror + PTTpos_mirror3;
            PTT_mirror = mapSegments(PTTpos_mirror);
            DM1.setIrisAO(PTT_mirror);
            
        else %Turn on correction again
            [ ~, PTTpos_mirror1] = IrisAOgetPTT_new( dOTF, pixelshift, lambda, [22,23,24], 23, calibration_filename_23);
            [ dOTF2, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF_new( DM1, 35, PTTpos_mirror, Noise_Parameters, F, A, ABER, TURB );
            dOTF2 = conj(dOTF2);
            [ PTT_mirror2, PTTpos_mirror2] = IrisAOgetPTT_new( dOTF2, pixelshift, lambda, [34,35,36], 35, calibration_filename_35);
            PTTpos_mirror3 = PTTpos_mirror1;
            PTTpos_mirror3([22,23,24],:) = PTTpos_mirror2([22,23,24],:);
            PTTpos_mirror3 = -PTTpos_mirror3;
            PTTpos_mirror3(:,1) = PTTpos_mirror3(:,1) * 0.654008264745021162;
            PTTpos_mirror = PTTpos_mirror + PTTpos_mirror3;
            PTT_mirror = mapSegments(PTTpos_mirror);
            DM1.setIrisAO(PTT_mirror);
        end
    end
    
    %store the commands sent to the DM
    DMCOMMANDSCUBE(:,:,nn) = PTT_mirror;
    
    %Get the Residual Field
    F.touch;
    F.planewave * ABER * TURB * A * DM1;

    %Compute the RMS Wavefront Error
    W = angle(F.grid).*hexmask;
%     W_sqar = W.^2;
%     W_sqar_mean = mean(W_sqar(abs(W_sqar)>0));
%     W_mean = mean(W(abs(W)>0));
%     W_mean_sqar = W_mean^2;
    %     WFE(nn) = sqrt(W_sqar_mean - W_mean_sqar);
    WFE(nn) = var(W(abs(W)>0));
    
    %Compute the Marechal Strehl Approximation
    strehl(nn) = exp(-WFE(nn));
    
    
    %Compute the Corrected PSF
    PSF_cor = F.mkPSF(FOV,PLATE_SCALE);
    maxPSF_cor = max(max(PSF_cor));
    if Noise_Parameters{5} == true
        PSF_cor = addNoise(PSF_cor,Noise_Parameters{6},Noise_Parameters{2},Noise_Parameters{3},Noise_Parameters{4});
        PSF_cor = abs(PSF_cor);
    end
    
    maxPSF_cor_noise = max(max(PSF_cor));

    %Compute the approximate strehls in the old way
%     strehl_cor(nn) = PSF_cor(ceil(length(PSF_cor)/2),ceil(length(PSF_cor)/2)) / PSF_difflim(ceil(length(PSF_difflim)/2),ceil(length(PSF_difflim)/2));
%     strehl_uncorr(nn) = PSF_aberrated(ceil(length(PSF_aberrated)/2),ceil(length(PSF_aberrated)/2)) / PSF_difflim(ceil(length(PSF_difflim)/2),ceil(length(PSF_difflim)/2));
%     strehl_uncorr(nn) = PSF_aberratedmax / PSF_difflimmax;
    strehl_notip(nn) = maxPSF_cor / PSF_difflimmax;
    loopnum(nn) = nn;
    
    
    
    
%     [DM1x,DM1y] = DM1.coords;
%     
    %     figure(1);
    %     subplot(2,4,1)
    %     imagesc(plotx,ploty,log10(PSF_difflim / PSF_difflimmax),[-4,0]);
    %     sqar;
    %     title('Diffraction Limited PSF');
    %     subplot(2,4,2)
    %     imagesc(plotx,ploty,log10(PSF_aberrated / PSF_aberratedmax),[-4,0])
    %     sqar;
    %     title('Aberrated PSF')
    %     subplot(2,4,5);
    %     imagesc(plotx,ploty,log10(PSF_cor / maxPSF_cor),[-4,0]);
    %     sqar;
    %     title(sprintf('Corrected PSF, loop #%d',nn));
    %     subplot(2,4,6);
    %     strehl(nn) = PSF_cor(401,401) / PSF_difflim(401,401);
    %     strehl_uncorr(nn) = PSF_aberrated(401,401) / PSF_difflim(401,401);
    %     strehl_notip(nn) = maxPSF_cor / PSF_difflimmax;
    %     loopnum(nn) = nn;
    %     plot(loopnum,strehl,'-r');
    %     hold on
    %     plot(loopnum,strehl_uncorr,'-b');
    %     plot(loopnum,strehl_notip,'-g');
    %     hold off
    %     xlabel('Loop Iteration');
    %     ylabel('Strehl Ratio');
    %     legend('Strehl for Corrected PSF','Strehl for Uncorrected/Aberrated PSF','Tip/Tilt Independent Strehl','Location','Best');
    %     xlim([0,numiterations]);
    %     ylim([0,1]);
    %     title('Strehl Ratio');
    %     %     drawnow;
    %     %
    %     %
    %     %     figure(2)
    %     subplot(2,4,3);
    %     imagesc(DM1x,DM1y,angle(DM1.grid));
    %     colorbar;
    %     axis xy;
    %     sqar;
    %     bigtitle(sprintf('Phase of DM with Correction Applied\n'),10);
    %     subplot(2,4,4);
    %     F.planewave * TURB * ABER * DM1;
    %     F.show;
    %     axis xy;
    %     sqar;
    %     colorbar;
    %     bigtitle('Field at DM',10);
    %     subplot(2,4,7)
    %     plotComplex(dOTF,3);
    %     axis off;
    %     axis xy;
    %     sqar;
    %     colorbar;
    %     bigtitle(sprintf('dOTF and Segment Center Locations\n'),10);
    %     hold on
    %     for n = 1:37
    %         if n ~= 23
    %             plot(pixel_seg_map{n}(2),pixel_seg_map{n}(1),'r*');
    %         else
    %         end
    %     end
    %     hold off
    %     subplot(2,4,8)
    %     if InjectKolm == true
    %         TURB.show;
    %         sqar;
    %         bigtitle(sprintf('Turbulence Profile at loop %d \n',nn),10);
    %         xlim([-2e-3,2e-3]);
    %         ylim([-2e-3,2e-3]);
    %     else
    %         ABER.show;
    %         bigtitle(sprintf('Injected Aberration at loop %d \n',nn),10);
    %         xlim([-2e-3,2e-3]);
    %         ylim([-2e-3,2e-3]);
    %         sqar;
    %     end
    
    
    subplot(2,2,1);
    imagesc(plotx,ploty,log10((PSF_cor / maxPSF_cor_noise)),[-4,0]);
    %     imagesc(plotx,ploty,PSF_cor);
    colormap(gray);
    sqar;
    colorbar;
    bigtitle(sprintf('PSF, loop #%d',nn),12);
    
    subplot(2,2,2);
    %     plot(loopnum,strehl_cor,'-r');
    %     hold on
    %     plot(loopnum,strehl_uncorr,'-b');
    %     plot(loopnum,strehl_notip,'-g');
    %     hold off
    %     xlabel('Loop Iteration');
    %     ylabel('Strehl Ratio');
    %     legend('Strehl for Corrected PSF','Strehl for Uncorrected/Aberrated PSF','Tip/Tilt Independent Strehl','Location','Best');
    %     xlim([0,numiterations]);
    %     ylim([0,1]);
    %     title('Strehl Ratio');
    
    hold on
    plot(loopnum,strehl,'-b');
    plot(loopnum,strehl_notip,'--r');
    plot(loopnum,strehl_uncorr,'--k');
    hold off;
    xlim([0,numiterations]);
    ylim([0,1]);
    bigtitle('Strehl Ratio',10);
    legend('Marechal Strehl Ratio of Correction Signal','Maximum Intensity Strehl Ratio of Correction Signal','Marechal Strehl Ratio of Injected Aberration Signal','Location','SouthOutside');
    
    
    subplot(2,2,3);
    F.planewave * TURB * ABER * DM1;
    F.show;
    axis xy;
    sqar;
    colorbar;
    bigtitle('Residual Field',12);
    
    subplot(2,2,4)
    plotComplex(dOTF,6);
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
    
    drawnow;
    MOVIE(:,nn) = getframe(moviefig,winsize);
    fprintf('Loop #%d Approximate Strehl: %0.6f \n',nn,strehl(nn));
    nn = nn + 1;
end


mkdir(filename);
current_dir = pwd;
cd(filename);

% fprintf('Saving OTFaCUBE\n');
% save(filename_OTFa,'OTFaCUBE','-v7.3');
% 
% fprintf('Saving OTFbCUBE\n');
% save(filename_OTFb,'OTFbCUBE','-v7.3');
% 
fprintf('Saving DM Positions\n');
save('PTT_commands_for_each_loop.mat','DMCOMMANDSCUBE','-v7.3');

fprintf('Saving Movie Variable\n');
save('Movie_and_Noise_Parameters.mat', 'MOVIE', 'winsize', 'Noise_Parameters','-v7.3')

fid = fopen('Simulation_Output_Data.txt','w+');
fprintf(fid,'Output Data from %s Run \r\n\n',filename);
fprintf(fid,'Pixel Shift Used: [%d %d] \r\n\n',pixelshift);

if Noise_Parameters{5} == true
    fprintf(fid,'Noise Parameters \r\n');
    fprintf(fid,'Number of Exposures per PSF: %d \r\n',Noise_Parameters{1});
    fprintf(fid,'Number of Photons per Exposure: %d \r\n',Noise_Parameters{6});
    fprintf(fid,'Shot Noise Included? %d \r\n',Noise_Parameters{2});
    fprintf(fid,'Standard Deviation of Read Noise: %0.2f \r\n',Noise_Parameters{3});
    fprintf(fid,'Dark Current Scale: %0.2f \r\n',Noise_Parameters{4});
end

if InjectAb == true
    fprintf(fid,'Zernike Aberrations Used: n, m, coefficient: \r\n');
    fprintf(fid,'%d \t',n_zern);
    fprintf(fid,'\r\n');
    fprintf(fid, '%d \t',m);
    fprintf(fid,'\r\n');
    fprintf(fid,'%0.4f \t',coeffs);
end
if InjectKolm == true
    fprintf(fid,'\r\n\n');
    fprintf(fid,'Using Kolmogorov Phase Screen\r\n');
    fprintf(fid,'r0 = %0.5f metes \r\n',r0);
    fprintf(fid,'Wind = [%d %d]',Wind(1),Wind(2));
end
fprintf(fid,'\r\nLoop Data \r\n');
fprintf(fid,'Correction Signal Data \r\n');
for n = 1:length(loopnum)
    fprintf(fid,'Loop %d: %0.5f Strehl, %0.5f WFE \r\n',loopnum(n),strehl(n),WFE(n));
end
fprintf(fid,'Uncorrected Signal Data \r\n');
for n = 1:length(loopnum)
    fprintf(fid,'Loop %d: %0.5f Strehl, %0.5f WFE \r\n',loopnum(n),strehl_uncorr(n),WFE_uncor(n));
end
fclose(fid);
cd(current_dir);

% movie(figure(1),MOVIE,1,10,winsize)
% movie2avi(MOVIE,filename_movie_avi,'compression','None','fps',7)




% Zernike_Number = [2];
% Zernike_Coefficient_waves = 4;
% PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);
% PTT = mapSegments(PTTpos);
% pistonlist_zern = PTT(:,1);
% pistonlist = PTT_mirror(:,1);
% pistonlist = horzcat(pistonlist,pistonlist_zern);
% tiplist = horzcat(PTT_mirror(:,2),PTT(:,2));
% tiltlist = horzcat(PTT_mirror(:,3),PTT(:,3));

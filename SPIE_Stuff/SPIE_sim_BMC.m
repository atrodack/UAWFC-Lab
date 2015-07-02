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
D = 7e-3; % 7mm
% secondary = 0.3*D; % 30% of Pupil Diameter
% spider = 0.02*D; % 2% of Pupil Diameter



%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 9; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
numiterations = 100;
gain = 0.7; %gain for AO Corrections
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
IrisAO_on = false; %turns on/off IrisAO Mirror (if false, DM1 variable set to 1)
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed
Scalloped_Field = true; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.

% BMC Flag
BMC_on = true; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

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



% make a field object
display('Making Fields');
new_spacing = DM2.spacing;
F = AOField(DM2);
F.FFTSize = (fftsize);
F.spacing(new_spacing);
F.lambda = lambda;
F.FOV = FOV;
F.PLATE_SCALE = PLATE_SCALE;
F.FoV = FoV_withIrisAO;
F.name = 'IrisAO Field 1';

% make diffraction limited PSF
display('Making Diffraction Limited PSF');
F.planewave * A * DM2;
[PSF_difflim,plotx,ploty] = F.mkPSF(FOV,PLATE_SCALE);
PSF_difflimmax = max(max(PSF_difflim));
F.touch;

%% Control Loop
nn = 1;

while(nn < numiterations)
    
    if InjectKolm == true
        TURB.grid(circshift(TURB.grid,Wind));
    end
    
    F.planewave * A * ABER * TURB * DM2;
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;
    


    
    %% Correction
    if nn > 10  %suffer the seeing limit for a bit
        
    else
        
    end
    
    F.touch;
    F.planewave * A * ABER * TURB * DM2;
    % figure(2);
    % F.show;
    
    PSF_cor = F.mkPSF(FOV,PLATE_SCALE);
    maxPSF_cor = max(max(PSF_cor));
    
    
    [DM2x,DM2y] = DM2.coords;
    
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
%     xlim([0,100]);
%     ylim([0,1]);
%     title('Strehl Ratio');
% %     drawnow;
% %     
% %     
% %     figure(2)
%     subplot(2,4,3);
%     imagesc(DM1x,DM1y,angle(DM2.grid));
%     colorbar;
%     axis xy;
%     sqar;
%     bigtitle(sprintf('Phase of DM with Correction Applied\n'),10);
%     subplot(2,4,4);
%     F.planewave * TURB * ABER * DM2;
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
%     
%     
%     drawnow;
%     
%     fprintf('Approximate Strehl: %0.5f \n',strehl(nn));
    nn = nn + 1;
end


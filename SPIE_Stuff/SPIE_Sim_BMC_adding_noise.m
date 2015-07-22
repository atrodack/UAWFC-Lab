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
% secondary  0.3*D; % 30% of Pupil Diameter
secondary = 0;
% spider = 0.02*D; % 2% of Pupil Diameter
spider = 0;


%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
nzerns = 4; %number of zernikes to inject (if InjectAb and InjectRandAb are both true)
% numiterations = 2;
gain = -0.6; %gain for AO Corrections
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
Scalloped_Field = false; %turns on/off returning an AOField Object that encodes the actual surface shape of the segments.

% BMC Flag
BMC_on = true; %turns on/off BMC Mirror (if false, DM2 variable is set to 1)

% Aberration Flags
InjectAb = true; %Injects nzerns Zernike Terms
InjectRandAb = false; %if InjectAB is true, picks Zernikes "Randomly"
InjectKnownAb = true; %if InjectAB is true, picks provided Zernikes

InjectKolm = true;
InjectKnownKolm = true;

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


% Compute the Number of Photons
Quantum_Efficiency = 0.5;
Bandpass = 0.1; %in microns
Exposure_Time = 0.200; %in seconds (approximate current testbed camera time)
Band_Flux = AOField.RBANDF; % in ph*um^-1*m^-2*s^-1
Star_Visual_Mag = 7;
D_Telescope = 6.5; %in meters (MMT)


N0 = Quantum_Efficiency * Bandpass * Exposure_Time * Band_Flux * (2.512^(-Star_Visual_Mag)) * ((pi*D_Telescope^2)/4);


if UseNoise == true
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 1;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 10;
    Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = N0;
else
    Noise_Parameters = cell(5,1);
    %     Noise_Parameters{1} = 1;
    %     Noise_Parameters{2} = false;
    %     Noise_Parameters{3} = 0;
    %     Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = N0;
end

% Diameter = 6.5;
% Area = pi * (Diameter/2)^2
% 10^10 * Area
% Q = 10^10 * Area * 0.2 * 0.33 * (1/1000)
% Zero_mag = Q;
% fith_mag = Zero_mag /100;


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
                n_zern = [1,4,4,3];
%         n_zern = [1];
        %                 m = [-2,0,2,-1,3];
                m = [1,-4,-2,1];
%         m = [1];
        
        %  coeffs = 1 * randn(1,length(n_zern));
        %  coeffs = [0.2441,-0.0886884,2.75*-0.0980274,-0.05,0.12];
        %  coeffs = 0.25*randn(1,length(n_zern));
        coeffs = [0.4575 	0.4649 	-0.3424 	0.4706 	];
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
        if InjectKnownKolm == false
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
        else
            load KnownTURB1.mat
            Wind = [-2,2];
        end
        r0 = TURB.r0;
        D_fit = D;
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
% if Noise_Parameters{5} == true
%     s1 = sprintf('%d_photons',Noise_Parameters{6});
%     s3 = sprintf('%0.2f_ReadNoise',Noise_Parameters{3});
%     s4 = sprintf('%d_images_per_PSF',Noise_Parameters{1});
% else
%     s1 = 'No_Noise_Used';
%     s3 = '_';
%     s4 = '_';
% end
% if InjectAb == true
%     s2 = sprintf('Zernike');
%     if InjectKolm == true
%         s2 = sprintf('Zernike_with_Kolmogorov_added');
%     end
% else
%     if InjectKolm == true
%         s2 = sprintf('Kolmogorov');
%     else
%         s2 = sprintf('No_Aberration');
%     end
% end


dt = datestr(now,'mm_dd_HH_MM_SS');
filename = sprintf('BMC_Deconv_Test_YesD_%s',dt) %add date and time to msec in order to avoid overwriting files
filename_movie_avi = sprintf('%s.avi',filename);


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

%set actuators that are off equal to empty
for n = 1:DM2.nActs
    if DM2.actuators(n,5) == 0
        calibrated_BMC_act_locations{n} = [];
    end
end
%removie the off actuators from the dOTF "reconstructor" list
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


Ppos = zeros(DM2.nActs,1);
Ppos2 = zeros(length(onAct_locations),1);

%% Control Loop
nn = 1;
no_actuators = [1,32,993,1024];
numiterations = 54;%60;
correction_start = 6;%6;

centerpoint = [746,854];
centerpoint2 = [1304,1196-9];
radius = 550;
radius2 = 350;
[X,Y] = meshgrid(1:fftsize);
R  = sqrt((X-centerpoint(2)).^2 + (Y - centerpoint(1)).^2);
SNRmask = double(R<=radius);
R2 = sqrt((X-centerpoint2(2)).^2 + (Y - centerpoint2(1)).^2);
CONJmask = ~double(R2<=radius2);
SNRmask = SNRmask .* CONJmask;

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
DMCOMMANDSCUBE = zeros(DM2.nActs,1);

strehl = zeros(1,1);
strehl_notip = zeros(1,1);
strehl_uncorr = zeros(1,1);
load dOTF_act_698_mask.mat;
load circmask.mat;

fprintf('Starting the loop\n\n');
DECONVOLVE = true;


DECONVOLVE = true;


while(nn <= numiterations)
    
    loopnum(nn) = nn;
    
    
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
    Ppos = DM2.actuators(:,3);
    DM2.flatten;
    F.planewave * TURB *ABER * A * DM2;
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;

    W = angle(F.grid).*circmask;
    WFE_uncor(nn) = var(W(abs(W)>0));
    strehl_uncorr(nn) = exp(-WFE_uncor(nn));
    
    PSF_aberrated = F.mkPSF(FOV,PLATE_SCALE);
    PSF_aberratedmax = max(max(PSF_aberrated));
    F.touch;
    
    
    [ dOTF, PSF1, PSF2, OTF1, OTF2, dOTFD ] = BMCcomputedOTF( DM2, 698, Ppos, Noise_Parameters, F, A, ABER, TURB, DECONVOLVE );

    dOTF = -1i * conj(dOTF);
    dOTFD = -1i * conj(dOTFD);
    
    if DECONVOLVE == true
        dOTF2 = dOTF;
        dOTF = dOTFD;
        clear dOTFD;
    end
    
    %Store the dOTF and PSF1 data for access later
%     PSFaCUBE(:,:,nn) = PSF1;
%     PSFbCUBE(:,:,nn) = PSF2;
%     OTFaCUBE(:,:,nn) = OTF1;
%     OTFbCUBE(:,:,nn) = OTF2;
%     dOCUBE(:,:,nn) = dOTF;
    azi_avg = azimuthal_average(fftshift(circshift(SNRmask .* abs(dOTF).^2,1-centerpoint)));
%     plot(log10(azi_avg/max(azi_avg)));
    S = mean((azi_avg(1:radius2)));
    N = mean((azi_avg(radius2:radius)));
    SNR(nn) = S/N;
    
    
    %% Correction
    if nn < correction_start  %suffer the seeing limit for a bit
        %Flatten the DM for the seeing limit
        DM2.flatten;
        
    elseif nn >= correction_start && nn < 3*correction_start %close the loop
        %Get Segment Piston from dOTF
        [ Ppos2, masked_dOTF ] = BMCgetP( dOTF, DM2,lambda, onAct_locations,'gold');
        
        %Bump the Actuator Pistons
        DM2.bumpOnActs(gain * (Ppos2));
        %Figure out the Slaves
        if nn == correction_start
            pistonlists = zeros(length(DM2.OnActs),2);
            pistonlists(:,1) = linspace(1,length(DM2.OnActs),length(DM2.OnActs));
            pistonlists(:,2) = Ppos2;
            slaveActs_dOTF_space = pistonlists(pistonlists(:,2)==0);
            slaveActs = DM2.OnActs(slaveActs_dOTF_space);
        end
        
        pistonlist(:,1) = DM2.actuators(:,3);
        Ppos2 = setOverlapActs(pistonlist,slaveActs);
        DM2.flatten;
        DM2.setActs(Ppos2);
        DM2.clip(STROKE);

        
        %Update the Mirror
        DM2.touch;
        DM2.render;
%         dOTF = masked_dOTF;
        
    elseif nn >= 3*correction_start && nn < 6*correction_start
        if InjectAb == true && InjectKolm == true %Simulate adding in some hair spray
            TURB = PHASESCREEN.copy;
            [ Ppos2, masked_dOTF ] = BMCgetP( dOTF, DM2, lambda, onAct_locations,'gold');
            DM2.bumpOnActs(gain * (Ppos2));
            pistonlist(:,1) = DM2.actuators(:,3);
            Ppos2 = setOverlapActs(pistonlist,slaveActs);
            DM2.flatten;
            DM2.setActs(Ppos2);
            DM2.clip(STROKE);
            DM2.touch;
            DM2.render;
%             dOTF = masked_dOTF;
                    
        else %Turn off correction for a bit
            %Flatten the DM
            DM2.flatten;
        end

%         
    elseif nn >= 6*correction_start
        if InjectAb == true && InjectKolm == true %Simulate taking the hair spray out
            TURB = 1;
            [ Ppos2, masked_dOTF ] = BMCgetP( dOTF, DM2, lambda, onAct_locations,'gold');
            DM2.bumpOnActs(gain * (Ppos2));
            pistonlist(:,1) = DM2.actuators(:,3);
            Ppos2 = setOverlapActs(pistonlist,slaveActs);
            DM2.flatten;
            DM2.setActs(Ppos2);
            DM2.clip(STROKE);
            DM2.touch;
            DM2.render;
%             dOTF = masked_dOTF;

        else %Turn on correction again
            [ Ppos2, masked_dOTF ] = BMCgetP( dOTF, DM2, lambda, onAct_locations,'gold');
            DM2.bumpOnActs(gain * (Ppos2));
            pistonlist(:,1) = DM2.actuators(:,3);
            Ppos2 = setOverlapActs(pistonlist,slaveActs);
            DM2.flatten;
            DM2.setActs(Ppos2);
            DM2.clip(STROKE);
            DM2.touch;
            DM2.render;
%             dOTF = masked_dOTF;

        end
    end
    
    %store the commands sent to the DM
    DMCOMMANDSCUBE(:,nn) = DM2.actuators(:,3);
    
    %Get the Residual Field
    F.touch;
    F.planewave * ABER * TURB * A * DM2;

    %Compute the RMS Wavefront Error
    W = angle(F.grid).*circmask;
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

    strehl_notip(nn) = maxPSF_cor / PSF_difflimmax;
    
    
    
    
    
    subplot(2,3,1);
    imagesc(plotx,ploty,log10((PSF_cor / maxPSF_cor_noise)),[-4,0]);
    %     imagesc(plotx,ploty,PSF_cor);
    colormap(gray);
    sqar;
    axis xy;
    axis off;
    colorbar;
    bigtitle(sprintf('PSF, loop #%d',nn),12);
    
    subplot(2,3,4);
    hold on
    plot(loopnum,strehl,'-b');
    plot(loopnum,strehl_notip,'--r');
    plot(loopnum,strehl_uncorr,'--k');
    hold off;
    xlim([0,numiterations]);
    ylim([0,1]);
    xlabel('Loop Iteration');
    ylabel('Strehl Value');
    bigtitle('Strehl Ratio',10);
    legend('Marechal Strehl Ratio of Correction Signal','Maximum Intensity Strehl Ratio of Correction Signal','Marechal Strehl Ratio of Injected Aberration Signal','Location','SouthOutside');
    
    
    subplot(2,3,3);
    F.planewave * TURB * ABER * A * DM2;
    F.show;
    axis xy;
    sqar;
    colorbar off;
    bigtitle('Residual Field',12);
    
    subplot(2,3,2)
    plotComplex(dOTF,6);
    axis off;
    axis xy;
    sqar;
    colorbar off;
    if DECONVOLVE == false
        bigtitle(sprintf('dOTF before Correction is Applied\n'),10);
    else
        bigtitle(sprintf('Deconvolved dOTF before Correction is Applied\n'),10);
    end
    
    %     hold on
    %     for n = 1:length(DM2.OnActs)
    %         plot(onAct_locations{n}(1),onAct_locations{n}(2),'g.');
    %     end
    %     hold off
    if DECONVOLVE == true
        subplot(2,3,5)
        plotComplex(dOTF2,6);
        axis off;
        axis xy;
        sqar;
        colorbar off;
        bigtitle(sprintf('Not Deconvolved dOTF before Correction is Applied\n'),10);
    end
    
    
    subplot(2,3,6)
    A_C = A.copy;
    A_C * DM2;
    A_C.show;
    colormap(jet);
    bigtitle(sprintf('DM Shape\n'));
    axis off
    

    drawnow;
    MOVIE(:,nn) = getframe(moviefig,winsize);
    fprintf('Loop #%d Approximate Strehl: %0.6f\t\tApproximate dOTF SNR: %0.4f \n',nn,strehl(nn),SNR(nn));
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
save('DM_commands_for_each_loop.mat','DMCOMMANDSCUBE','-v7.3');

fprintf('Saving Movie Variable\n');
save('Movie_and_Noise_Parameters.mat', 'MOVIE', 'winsize', 'Noise_Parameters','-v7.3')

fid = fopen('Simulation_Output_Data.txt','w+');
fprintf(fid,'Output Data from %s \r\n\n',filename);

fprintf(fid,'Quantum Efficiency: %0.2f \r\n',Quantum_Efficiency);
fprintf(fid,'Bandpass: %0.4f in microns \r\n',Bandpass);
fprintf(fid,'Exposure Time: %0.3f in sec \r\n',Exposure_Time);
fprintf(fid,'Band Flux: %d in ph * s^-1 * m^-2 * um^-1 \r\n',Band_Flux);
fprintf(fid,'Apparent Star Magnitude: %0.2f \r\n',Star_Visual_Mag);
fprintf(fid,'Telescope Diameter: %0.3f meters \r\n',D_Telescope);
fprintf(fid,'\r\n');
if Noise_Parameters{5} == true
    fprintf(fid,'Noise Parameters \r\n');
    fprintf(fid,'Number of Exposures per PSF: %d \r\n',Noise_Parameters{1});
    fprintf(fid,'Number of Photons per Exposure: %d \r\n',Noise_Parameters{6});
    fprintf(fid,'Shot Noise Included? %d \r\n',Noise_Parameters{2});
    fprintf(fid,'Standard Deviation of Read Noise: %0.2f \r\n',Noise_Parameters{3});
    fprintf(fid,'Dark Current Scale: %0.2f \r\n',Noise_Parameters{4});
end
fprintf(fid,'\r\n');
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
    fprintf(fid,'Loop %d: %0.5f Strehl, %0.5f WFE, %0.5f SNR \r\n',loopnum(n),strehl(n),WFE(n),SNR(n));
end
fprintf(fid,'\r\nUncorrected Signal Data \r\n');
for n = 1:length(loopnum)
    fprintf(fid,'Loop %d: %0.5f Strehl, %0.5f WFE \r\n',loopnum(n),strehl_uncorr(n),WFE_uncor(n));
end
fclose(fid);
cd(current_dir);

% movie(figure(1),MOVIE,1,10,winsize)
% movie2avi(MOVIE,filename_movie_avi,'compression','None','fps',7)


















































































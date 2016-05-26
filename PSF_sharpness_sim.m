%% Start Clean
clear all;
clc;
clf;
% close all;


% Set Wavelength
% lambda = AOField.RBAND; % Red light.
lambda = AOField.HeNe_Laser;

% Compute Wavenumber
k = (2*pi) / lambda;

% Pupil Specs
D = 7e-3; % 7mm
secondary = 0;
secondary = 0.3*D; % 30% of Pupil Diameter
spider = 0;
spider = 0.02*D; % 2% of Pupil Diameter

%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 7.5*SPACING;  % for antialiasing.
fftsize = 2^11;

%% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 10*THld; % arcsecs
PLATE_SCALE = THld/10;
FoV_withoutIrisAO = 10.5e-3;


% Compute the Number of Photons
Quantum_Efficiency = 0.5;
Bandpass = 0.1; %in microns
Exposure_Time = 0.200; %in seconds (approximate current testbed camera time)
Band_Flux = AOField.RBANDF; % in ph*um^-1*m^-2*s^-1
Star_Visual_Mag = 7;
D_Telescope = 6.5; %in meters (MMT)


N0 = Quantum_Efficiency * Bandpass * Exposure_Time * Band_Flux * (2.512^(-Star_Visual_Mag)) * ((pi*D_Telescope^2)/4);


%% Make Testbed Elements
RunSIM=true;
RunTESTBED=false;
system_verbose = false;
IrisAO_on = false;
BMC_on = true;
maketestbedelements;
DM = DM2.copy;
DM.name = 'Kilo DM';
clearvars -except DM lambda k D SPACING aa fftsize THld FOV PLATE_SCALE UseNoise Noise_Parameters N0 A A_BMC;
clc;

% Noise Flags
UseNoise = false;

if UseNoise == true
    fprintf('******* Noise is being included! *******\n\n');
    Noise_Parameters = cell(5,1);
    Noise_Parameters{1} = 1;
    Noise_Parameters{2} = true;
    Noise_Parameters{3} = 10;
    Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = N0;
else
    Noise_Parameters = cell(5,1);
    %         Noise_Parameters{1} = 1;
    %     Noise_Parameters{2} = false;
    %     Noise_Parameters{3} = 0;
    %     Noise_Parameters{4} = 0;
    Noise_Parameters{5} = UseNoise;
    Noise_Parameters{6} = N0;
end


%% Plot Pupil Mask and DM Acts/Regions
system_verbose = false;
if system_verbose == true;
    figure(1);
    subplot(1,2,1);
    A.show; sqar; colormap(gray(256));
    subplot(1,2,2);
    DM.show; sqar; colormap(gray(256));
    pause(1);
    DM.plotActuators;
    pause(1);
    DM.plotRegions;
end
%% Make a Field
F = AOField(A);
F.spacing(SPACING);
F.FFTSize = fftsize;
% F.resize(fftsize);
F.lambda = lambda;
F.planewave;
F*A*DM;
F.show
[PSF_Difflim,THX,THY] = F.mkPSF(FOV,PLATE_SCALE);
% imagesc(THX,THY,log10(normalize(PSF)),[-5,0]);sqar; colormap(gray(256));
F.touch;
if UseNoise == true;
    [   PSF_Difflim,~] = average_noisy_images(PSF_Difflim,0,N0,Noise_Parameters);
end

Diffraction_limited_norm = norm(PSF_Difflim,2);
fprintf('Norm for Diffraction Limited PSF: %g\n\n',norm(PSF_Difflim,2));

%% Add Zernikes to test norm merit function
% Nmax = 3;
% ABER = AOScreen(A);
% nmode = 1;
% for n=1:Nmax
%     fprintf('\nOrder %d:',n);
%     %     amp = 1*lambda/4/n;
%     amp = 1*lambda;
%     for m=-n:2:n
%         weight = (n^2 + m^2).^(-5/6);
%         ABER.zero.addZernike(n,m,weight*amp,D);
%         DM.setActs(ABER);
%         DM.render;
%         F.planewave * A * DM;
%         [PSF,THX,THY] = F.mkPSF(FOV,PLATE_SCALE);
%         [Noisy_PSF,~] = average_noisy_images(PSF,0,N0,Noise_Parameters);
%
%         %         imagesc(THX,THY,log10(normalize(PSF)),[-5,0]);sqar; colormap(gray(256));
%         imagesc(THX,THY,Noisy_PSF); sqar; colormap(gray(256));
%         F.touch;
%         DM.touch;
%         drawnow;
%         norm_PSF = norm(Noisy_PSF,2);
%         merit_func = norm_PSF/Diffraction_limited_norm;
%         nmode = nmode + 1;
%         fprintf('%d  %g\t ',m,merit_func);
%     end
% end
%
% fprintf('\n');


%% Add Zernikes to simulate our DM shape;
DM_sag = AOScreen(A);
DM_sag.zero;
amp = [-0.75,-0.85,0.41,0.66,0.75,0.5]*lambda;
n = [2,2,2,3,3,4];
m = [0,-2,2,-1,1,0];

for ii = 1:3
    weight = (n(ii)^2 + m(ii)^2).^(-5/6);
    DM_sag.addZernike(n(ii),m(ii),weight*amp(ii),D);
end


DM_sag.show;
DM.setActs(DM_sag);
DM.render;
F.planewave * A * DM;
[PSF,THX,THY] = F.mkPSF(FOV,PLATE_SCALE);
if UseNoise == true
    [PSF,~] = average_noisy_images(PSF,0,N0,Noise_Parameters);
end

%         imagesc(THX,THY,log10(normalize(PSF)),[-5,0]);sqar; colormap(gray(256));
imagesc(THX,THY,PSF); sqar; colormap(gray(256));
F.touch;
DM.touch;
drawnow;
norm_DM_PSF = norm(PSF,2);
merit_func = norm_DM_PSF/Diffraction_limited_norm;
fprintf('Ratio of norm with sag to Diffraction limit: %g\n',merit_func);
norm_vec = zeros(1,9);
amps = zeros(1,7);

%% Scan Zernikes
Nmax = 3;
ABER = AOScreen(A);
CORRECTION = AOScreen(A);
% CORRECTION.grid(ones(length(CORRECTION.grid)));
CORRECTION.zero;
nmode = 1;
counter2 = 1;
norm_PSF_prev = 0;
norm_PSF = 1;
for n=2:Nmax
    for m=-n:2:n
        if m == -3
            fprintf('Not doing m = -3\n');
        elseif m == 3
            fprintf('Not doing m = 3\n');
        else
            
            counter = 1;
            amp = -lambda;
            for amp = -lambda:0.01*lambda:lambda
                %         while(norm_PSF_prev < norm_PSF)
                norm_PSF_prev = norm_PSF;
                %             amp = amp + 0.1*lambda;
                weight = (n^2 + m^2).^(-5/6);
                ABER.zero.addZernike(n,m,weight*amp,D);
                DM.setActs(ABER);
                DM.render;
                F.planewave * A * DM_sag * DM;
                [PSF_zern,THX,THY] = F.mkPSF(FOV,PLATE_SCALE);
                if UseNoise == true
                    [PSF_zern,~] = average_noisy_images(PSF_zern,0,N0,Noise_Parameters);
                end
                %             imagesc(THX,THY,Noisy_PSF_zern); sqar; colormap(gray(256));
                F.touch;
                DM.touch;
                %             drawnow;
                norm_PSF = norm(PSF_zern,2);
                norm_vec(1,counter) = norm_PSF;
                merit_func = norm_PSF/norm_DM_PSF;
                nmode = nmode + 1;
                counter = counter + 1;
                fprintf('. ');
                
            end
            for ii = 1:length(norm_vec)
                if isequal(norm_vec(ii),max(norm_vec))
                    amp = -lambda:0.01*lambda:lambda;
                    amps(1,counter2) = amp(ii);
                end
                amp = amps(counter2);
            end
            %         amps(1,counter2) = amp;
            CORRECTION.addZernike(n,m,weight*amp,D);
            counter2 = counter2 + 1;
            fprintf('\n');
            norm_PSF_prev = 0;
            norm_PSF = 1;
        end
    end
end


CORRECTION.grid(-CORRECTION.grid);
F.planewave * A * DM_sag * CORRECTION * DM;

[PSF_corrected,THX,THY] = F.mkPSF(FOV,PLATE_SCALE);
if UseNoise == true
    [PSF_corrected,~] = average_noisy_images(PSF_corrected,0,N0,Noise_Parameters);
end
imagesc(THX,THY,PSF_corrected); sqar; colormap(gray(256));
F.touch;
DM.touch;
drawnow;
norm_corrected_PSF = norm(PSF_corrected,2);
























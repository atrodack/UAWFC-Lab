function dOTF = TestbeddOTF(DM, pokeact, Ppos_in)

if nargin < 3
    DM.flatten;
    Ppos_in = DM.actuators(:,3);
end

%% Initial Setup

Ppos_flat = Ppos_in;
DM.setActs(Ppos_flat);
DM_Pistons = reshape(DM2.actuators(:,3),[32,32]);
DM_Pistons = single(DM_Pistons);

tempdir = pwd;
cd /home/alex/Desktop/Testbed_fits_files;
fitswrite(DM_Pistons,'DM_Pistons.fits');
img = fitsread('DM_Pistons.fits');

figure(5);
imagesc(img);
input('Press Enter');
close
cd(tempdir);

% DO STUFF TO SEND TO MIRROR


% TAKE PICTURES ON TESTBED
% call run_cam

% PSF == Final Picture from here

input('Press Enter to Move On');

%% Modified Setup
Ppos_poked = Ppos_flat;
Ppos_poked(pokeact) = Ppos_poked(pokeact) + (AOField.HeNe_Laser*10^6) / 4;

DM.setActs(Ppos_poked);
DM_Pistons_poked = reshape(DM.actuators(:,3),[32,32]);
DM_Pistons_poked = single(DM_Pistons_poked);

cd /home/alex/Desktop/Testbed_fits_files;
fitswrite(DM_Pistons_poked,'DM_Pistons_poked.fits');
img = fitsread('DM_Pistons_poked.fits');
figure(5);
imagesc(img);
input('Press Enter');
close
cd(tempdir);

% DO STUFF TO SEND TO MIRROR


% TAKE PICTURES ON TESTBED
%call run_cam

% PSF_poked == Final Picture from here


%% Edit the Pictures
%shift to center, 






otf1 = fftshift(fft2(fftshift(PSF)));
otf2 = fftshift(fft2(fftshift(PSF_poked)));

dOTF = otf1 - otf2;

end
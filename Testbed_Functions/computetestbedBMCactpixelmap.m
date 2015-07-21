function calibrated_BMC_act_locations = computetestbedBMCactpixelmap(DM,pokeact)


DM.flatten;
DM = dmsetrow(DM,1,1e-6);
Ppos_flat = DM.actuators(:,3);
DM = dmsetcol(DM,1,1e-6);
Ppos_flat = Ppos_flat + DM.actuators(:,3);
Ppos_flat(Ppos_flat~=0) = 1;

DM.setActs(Ppos_flat);
DM_Pistons = reshape(DM.actuators(:,3),[32,32]);
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

% PSF == Final Picture from here

input('Press Enter to Move On');

Ppos_poked = Ppos_flat;
Ppos_poked(pokeact) = (AOField.HeNe_Laser*10^6) / 4;

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


% PSF_poked == Final Picture from here


otf1 = fftshift(fft2(fftshift(PSF)));
otf2 = fftshift(fft2(fftshift(PSF_poked)));

dOTF = otf1 - otf2;


figure;
plotComplex(dOTF,5);
axis xy;
sqar;



fprintf('Pick the Bottom Left Point\n\n');
btleft = pickPoint;
hold on
plot(btleft(2),btleft(1),'r*');
hold off


fprintf('Pick the Top Left Point\n\n');
topleft = pickPoint;
hold on
plot(topleft(2),topleft(1),'r*');
hold off


fprintf('Pick the Bottom Right Point\n\n');
btright = pickPoint;
hold on
plot(btright(2),btright(1),'r*');
hold off

x = linspace(btleft(2),btright(2),32);
y = linspace(btleft(1),topleft(1),32);
[X,Y] = meshgrid(x,y);


DM.flatten;
% center_col = ;
DM = dmsetcol(DM,center_col,1);
Ppos_flat = DM.actuators(:,3);

DM.setActs(Ppos_flat);
DM_Pistons = reshape(DM.actuators(:,3),[32,32]);
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

% PSF == Final Picture from here

input('Press Enter to Move On');

Ppos_poked = Ppos_flat;
Ppos_poked(pokeact) = (AOField.HeNe_Laser*10^6) / 4;

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


% PSF_poked == Final Picture from here

otf1 = fftshift(fft2(fftshift(PSF)));
otf2 = fftshift(fft2(fftshift(PSF_poked)));

dOTF = otf1 - otf2;


clf;
plotComplex(dOTF,5);
axis xy;
sqar;


% hold on
% plot(ceil(X),ceil(Y),'r*')
% hold off

calibrated_BMC_act_locations = cell(DM.nActs,1);
counter = 1;
for n = 1:32
    for m = 1:32
        hold on
        calibrated_BMC_act_locations{counter} = [X(m,n),Y(m,n)];
        calibrated_BMC_act_locations{counter} = round(calibrated_BMC_act_locations{counter});
        plot(calibrated_BMC_act_locations{n}(1),calibrated_BMC_act_locations{n}(2),'*g');
        drawnow;
        hold off
        counter = counter + 1;
    end
end



end
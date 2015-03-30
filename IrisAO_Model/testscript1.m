%% Clear Workspace
clear all;
clc;
close all;

%% Set Initial Parameters and Flags
segpitch = 606e-6;
magnification = 2;

verbose = false;
Scalloped_field = false;
lambda = AOField.RBAND;

%% Make the Mirror
if Scalloped_field == true
    [DM,F_Scal] = makeIrisAODM(magnification,verbose,Scalloped_field);
    F_Scal.name = 'Scalloped Field';
else
    DM = makeIrisAODM(magnification,verbose,Scalloped_field);
end


clc;
[x,y] = DM.coords;

extentx = (abs(min(x)) + abs(max(x)));
extenty = (abs(min(y)) + abs(max(y)));
% DM.Offset = [-x(ceil(length(x)/2)), 0];
fprintf('X Width = %0.4f mm\n',extentx * 10^3);
fprintf('Y Width = %0.4f mm\n',extenty * 10^3);

%% Make an Aperture and a Field
A = AOSegment(DM);
D = (extentx);
dx = (magnification *segpitch) / 100;
PNECO = [0 0 D 1 1.5*dx 0 0 0 0 0];
A.pupils = PNECO;
A.make;

DM.trueUp;
A.centerOn(DM);

F = AOField(A);
F.name = 'Centered IrisAO DM';
F.FFTSize = [2048 2048];
F.lambda = lambda;

%%
if Scalloped_field == false
load System_Tip_Correction_PTTpos.mat;
PTTpos = PTTPositionArray;
% DM = IrisAOPTT(DM,1:37,PTTpos(:,1)*10^-6,PTTpos(:,2)*10^-3,PTTpos(:,3)*10^-3,false);
% DM = IrisAOPTT(DM,1:37,zeros(37,1),PTTpos(:,2)*10^-3,PTTpos(:,3)*10^-3,false);
DM = IrisAOPTT(DM,1:37,zeros(37,1),zeros(37,1),zeros(37,1),false);
DM.show;
DM = IrisAOPTT(DM,1:37,zeros(37,1),zeros(37,1),ones(37,1)*-10^-3,true);
figure; DM.show;
DM = IrisAOPTT(DM,1:37,zeros(37,1),ones(37,1) *-10^-3,zeros(37,1),true);
figure; DM.show;
DM.grid;
F.planewave * DM * A;
F.show
colormap(gray);

else
    F_Scal * DM * A;
    F_Scal.show;
end

%% Make a PSF
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 50*THld; % FOV for PSF computation
PLATE_SCALE = THld/5; % Pixel Size for PSF computation
[PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
PSFmax = max(PSF(:));

if verbose == true
    figure(2)
    colormap(gray);
end
imagesc(thx,thy,log10(PSF/PSFmax),[-4,0]);
% imagesc(thx,thy,PSF/PSFmax);
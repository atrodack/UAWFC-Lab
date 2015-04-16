clear all;
clc;
close all;

%% UAWFC Testbed Simulation
lambda = AOField.RBAND; % Red light.
D = 7e-3; % 7mm
secondary = 0.3*D; % 30% of Pupil Diameter

SPACING = 1e-5; % fine spacing
aa = SPACING;  % for antialiasing.
spider = 0.02*D; % 2% of Pupil Diameter

% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 300; % arcsecs
PLATE_SCALE = THld/15;

% Flags
InjectAb = true;


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
fprintf('Pupil Mask Constructed\n');
figure(1)
A.show;
colormap(gray);
drawnow;
pause(1);

[x,y] = A.coords;

%% IrisAO DM
segpitch = 606e-6;
magnification = 1;

%Flags
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed

Scalloped_Field = true; %turns on/off returning an AOField Object that 
%encodes the actual surface shape of the segments.

% Make the Mirror
if Scalloped_Field == true
    [IrisAO_DM,F_scal] = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
else
    IrisAO_DM = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
end
IrisAO_DM.lambdaRef = lambda;

% clc; %clears the command window of the text generated from adding segments to the AOAperture object
fprintf('IrisAO Mirror Constructed\n');


% Set a Random Mirror Shape
% PTTpos = horzcat(horzcat(randn(37,1)*10^-6,randn(37,1)*10^-3),randn(37,1)*10^-3);

% Use a Zernike on the mirror
Zernike_Number = [2,3];
Zernike_Coefficient_waves = randn(2,1);
PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);


%
load('IrisAO_SegMap.mat');
PTT = zeros(37,3);

for ii = 1:37
    mapped_segment = IrisAO_SegMap(ii);
    PTT(ii,1:3) = PTTpos(mapped_segment,:);
end

IrisAO_DM.PTT(PTT);
IrisAO_DM.touch;
IrisAO_DM.render;
figure(1);
IrisAO_DM.show;
colormap(gray);
drawnow;
pause(1);

%% BMC DM
%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
Max_Stroke = 1.5e-6;
Pitch = 300e-6;
radius_BMC = (4.65)*10^-3;

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
A_BMC.show;
drawnow;

%Make it an AODM
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
BMC_DM.show;
colormap(gray);
title('BMC Pupil');
BMC_DM.plotActuators;
drawnow;

pause(1);
% BMC_DM.plotRegions;
% pause(1);
% fprintf('BMC DM Constructed\n');




%% Inject Random Zernike Aberration to simulate system errors

if InjectAb == true
    ABER = AOScreen(A);
    nzerns = 2;
    n = sort((randi(4,1,nzerns)),'ascend');
    m = zeros(1,nzerns);
    for ii = 1:nzerns
        m_pos = -n(ii):2:n(ii);
        choice = randi(length(m_pos),1,1);
        m(ii) = m_pos(choice);
    end
    ABER.zero;
    for ii = 1:nzerns
        ABER.addZernike(n(ii),m(ii),randn(1)*lambda,D);
    end
end

%% Go Through the System
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


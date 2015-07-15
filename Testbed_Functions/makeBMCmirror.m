








%% Simulation Parameters
SPACING = 1e-5; % fine spacing
aa = 5*SPACING;  % for antialiasing.
fftsize = 2^11;






%% Boston MircroMachines DM
%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
STROKE = 1.5e-6;
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

%Make an AODM
DM = AODM(A_BMC);
[X,Y] = DM.COORDS;

% Create Actuator Locations
actuator_x = xmin:Pitch:xmax;
actuator_x = actuator_x - mean(actuator_x);
actuator_y = ymin:Pitch:ymax;
actuator_y = actuator_y - mean(actuator_y);
[ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
%     ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
%     ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
%     ACTUATOR_X(ACTUATOR_X==0) = [];
%     ACTUATOR_Y(ACTUATOR_Y==0) = [];
BMC_ACTS = [ACTUATOR_X(:),ACTUATOR_Y(:)];
nacts = length(BMC_ACTS(:,1));

% Add Actuators
DM.addActs(BMC_ACTS,1,A_BMC);
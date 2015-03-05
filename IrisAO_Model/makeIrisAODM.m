function [varargout] = makeIrisAODM(MAGNIFICATION, SegPitch, Scalloped_Field)
% [DM] = makeIrisAODM(MAGNIFICATION, ACTUAL_D, SegPitch)
% Makes a 37 Segment IrisAO DM. SegPitch sets the Segment "Diameter",
% Magnification scales this, Scalloped_Field also returns a Scalloped Field

if nargin == 0
    MAGNIFICATION = 1;
    SegPitch = MAGNIFICATION * 606e-6;
    Scalloped_Field = false;
elseif nargin == 2
    SegPitch = MAGNIFICATION * SegPitch;
    Scalloped_Field = false;
elseif nargin == 3
    SegPitch = MAGNIFICATION * SegPitch;
else
    error('Inputs are Incorrect');
end

dx = SegPitch/100;
GAP = SegPitch*0.02;

Seg = AOSegment(ceil((SegPitch+2*GAP)/dx*1.25));
Seg.spacing(dx);

SCALLOPING = 100e-9;

[x,y] = Seg.coords;
[X,Y] = Seg.COORDS;

Q = (SegPitch/2)*exp(2*pi*1i*(0:6)'/6);
V = [real(Q) imag(Q)];
dV = diff(V,1,1);

SEG = 1;

for n=1:6
    th = 60/180*pi*n;
    SEG=SEG.*smoothedge(SegPitch/2- (X*cos(th)-Y*sin(th)),dx);
    imagesc(x,y,SEG);sqar;axis xy;
    colormap(gray);
    drawnow;
    pause(0.2)
end

Seg.grid(SEG);
Seg.name = 'IrisAO DM Segment';

DM = AOAperture();
Seg.name = 'IrisAO DM';
DM.spacing(dx);

a2 = (SegPitch/2+GAP)*1.333*1.0;
a1 = a2*sqrt(3)/2;

ua = a1 * [ 0  1.7 ];
th1 = 60/180*pi;
c = cos(th1); s = sin(th1);
ub = ([c -s; s c] * ua')';
ub_ = ([c s; -s c] * ua')';

figure(2)
for nrow=-3:0
    NumCol = 7-abs(nrow);
    
    START = -3*ua - nrow*ub;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
        DM.show;
        drawnow;
    end
end

for nrow=1:3
    NumCol = 7-abs(nrow);
    START = -3*ua + nrow*ub_;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
        DM.show;
        drawnow;
    end
end


SEGCOORDS = zeros(37,2);
for n=1:37
    SEGCOORDS(n,:) = DM.segList{n}.Offset;
end
varargout{1} = DM;
%% Scalloping on Field
if Scalloped_Field == true;
    F = AOField(DM);
    F.name = 'Science Field';
    F.lambda = AOField.HeNe_Laser;
    F.FFTSize = 1024*2;
    
    DM.lambdaRef = F.lambda;
    
    DM.touch;
    
    PS = AOScreen(DM);
    [X,Y] = PS.COORDS;
    
    PS.zero;
    w = SegPitch/2.2;
    
    for n=1:length(DM.segList)
        LOC = DM.segList{n}.Offset;
        R = sqrt((X-LOC(2)+SegPitch).^2+(Y-LOC(1)).^2); % Kludge!
        PS.grid(PS.grid+exp(-(R/w).^2));
    end
    
    PS.grid(SCALLOPING*normalize(PS.grid-min(PS.grid_(:))));
    
    DM.trueUp.touch.make;
    F.planewave*DM.touch.make*PS;
    F.show;
    varargout{2} = F;
end
end

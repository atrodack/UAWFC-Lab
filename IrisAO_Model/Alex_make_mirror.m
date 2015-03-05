% Make the IrisAO Mirror
clear all;
clc;
close all;


NN = 2;
dx = .01/NN;
GAP = 0.005;

Seg = AOSegment(64*NN);
Seg.spacing(dx);

[x,y] = Seg.coords;
[X,Y] = Seg.COORDS;

Q = sqrt(3)*max(x)*exp(2*pi*1i*(0:6)'/6);
V = [real(Q) imag(Q)];
dV = diff(V,1,1);

SEG = 1;

for n=1:6
    SEG=SEG.*smoothedge(dV(n,1)*(X-V(n,1))+dV(n,2)*(Y-V(n,2)),.003);
    %     imagesc(x,y,SEG);sqar;axis xy;
    %     colormap(gray);
    %     drawnow;
end

Seg.grid(SEG);
Seg.name = 'IrisAO DM Segment';

DM = AOAperture();
Seg.name = 'IrisAO DM';
DM.spacing(dx);
BB = Seg.BBox
a2 = max(abs(BB(:,2))) + GAP
a1 = max(abs(BB(:,1))) + GAP

ua = a1 * [ 0  1.7 ]
th1 = 60/180*pi;
c = cos(th1); s = sin(th1);
ub = ([c -s; s c] * ua')'
ub_ = ([c s; -s c] * ua')'

for nrow=-3:0
    NumCol = 7-abs(nrow)
    
    START = -3*ua - nrow*ub;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
    end
end

for nrow=1:3
    NumCol = 7-abs(nrow)
    START = -3*ua + nrow*ub_;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
    end
end


SEGCOORDS = zeros(37,2);
for n=1:37
        SEGCOORDS(n,:) = DM.segList{n}.Offset;
end

DM.show;

fprintf('\n\n\n');




%% Put a Field on It
A = AOSegment(DM);
D = 1.75;
% PNECO = [0 0 3.4 1 0.05 0 0 0 0 0];
PNECO = [0 0 D 1 0.02 0 0 0 0 0];
A.pupils = PNECO;
A.make;

[x,y] = A.coords;
% [X,Y] = Seg.COORDS;

DM.trueUp;

F = AOField(DM);
F.name = 'Science Field';
F.lambda = 650e-9;
F.FFTSize = 2^13;

DM.lambdaRef = F.lambda;
 
%% PTT
    segList = 1:37;
    
    %Piston
    % pistonList = rand(37,1)*10^-6;
    pistonList = zeros(37,1);
%     pistonList = randn(37,1)*10^-6;
    
%     Tip
    tipList = zeros(37,1);
%     tipList = randn(37,1)*10^-3;
    
%     Tilt
    tiltList = zeros(37,1);
    tiltList(1) = 10e-3;
    tiltList(6) = 10e-3;

%     tiltList = randn(37,1)*10^-3;

    DM = IrisAOPTT(DM,segList,pistonList,tipList,tiltList);
    
%% Put a Field on It
F.planewave*DM*A;
F.show;
drawnow;

FOV=5
dFOV = FOV/1000
PSF = F.mkPSF(FOV,dFOV);
figure(2)
subplot(1,2,1);
imagesc(PSF);
colormap(gray);
axis off;
subplot(1,2,2);
imagesc(log10(PSF));
axis off
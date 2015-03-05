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

imagesc(DM.grid);sqar;

A = AOSegment(DM);
D = 3.2;
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
F.FFTSize = 1024*2;

DM.lambdaRef = F.lambda;
 
% DM.segList{10}.piston = 150e-9
% DM.segList{20}.tiptilt = [1 2]/2000/NN;

DM.touch;
F.planewave*DM*A;
F.show;

FOV=1
dFOV = FOV/100
PSF = F.mkPSF(FOV,dFOV);

SCALE = .2

PSF_ = 0;
FCUBE  = zeros([F.size 100]);
IGRAMS = zeros([F.size 100]);

% 
% for loop=1:1
%     loop
%     DM.trueUp;
%     for n=1:37
%         DM.segList{n}.tiptilt = randn(1,2)/206265*SCALE;
%     end
%     % fix the pistons so the segments are phased...
%     DM.touch.make;
%     
%     PISTONS = angle(DM.interpGrid(SEGCOORDS(:,2),SEGCOORDS(:,1)))*DM.lambdaRef/2/pi;
%     DM.bumpPistons(-PISTONS/2);
%         
%     
%     %bar(PISTONS*1e9);
%     %ylim([-1 1]*1000);
%     % drawnow;
%     
%     F.planewave*DM.touch.make*A;
% 
%     FCUBE(:,:,loop) = F.grid;
%     IGRAMS(:,:,loop) = F.interferometer(1);
%     
%     subplot(1,2,2);
%     imagesc(x,y,IGRAMS(:,:,loop),[0 4]);sqar;
%     
%     subplot(1,2,1);
%     plotCAmpl(mean(FCUBE(:,:,1:loop),3),1/2);
% 
%     %imagesc(x,y,F.interferometer(1),[0 2]);sqar;
%     drawnow;
%     PSF = F.mkPSF(FOV,dFOV);
%     
%     PSF_ = PSF_ + PSF;
%     
% end
% 
% % MTF = abs(computeOTF(PSF));
% % MTF = abs(computeOTF(PSF_));
% 
% % fits_write_image('PSF.fits',PSF);
% % fits_write_image('PSF.fits',PSF);
% % fits_write_image('PSF_.fits',normalize(PSF_));
% % fits_write_image('MTF.fits',MTF);
% % fits_write_image('Interferograms.fits',IGRAMS);
% 
% 
% imagesc(log10(normalize(PSF)),[-5 0]);sqar;colorbar;
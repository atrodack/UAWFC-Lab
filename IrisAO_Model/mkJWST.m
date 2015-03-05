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

JWST = AOAperture();
Seg.name = 'IrisAO DM';
JWST.spacing(dx);
BB = Seg.BBox
a2 = max(abs(BB(:,2))) + GAP
a1 = max(abs(BB(:,1))) + GAP

ua = a1 * [ 0  1.7 ]
th1 = 60/180*pi;
c = cos(th1); s = sin(th1);
ub = ([c -s; s c] * ua')'
ub_ = ([c s; -s c] * ua')'

for nrow=-2:0
    NumCol = 5-abs(nrow)
    
    START = -2*ua - nrow*ub;
    
    for ncol=1:NumCol
        JWST.addSegment(Seg,START+ncol*ua);
    end
end

for nrow=1:2
    NumCol = 5-abs(nrow)
    START = -2*ua + nrow*ub_;
    
    for ncol=1:NumCol
        JWST.addSegment(Seg,START+ncol*ua);
    end
end


SEGCOORDS = zeros(length(JWST.segList),2);
for n=1:length(JWST.segList)
        SEGCOORDS(n,:) = JWST.segList{n}.Offset;
end

JWST.removeSegment(10);

imagesc(JWST.grid);sqar;

A = AOSegment(JWST);
PNIRCAM = [0 0 2*1.5 1 0.01 0 0 0 0 0];
A.pupils = PNIRCAM;
A.make.show;

[x,y] = A.coords;
% [X,Y] = Seg.COORDS;

JWST.trueUp;

F = AOField(JWST);
F.name = 'Science Field';
F.lambda = 650e-9;
F.FFTSize = 1024*2;

JWST.lambdaRef = F.lambda;
 
% DM.segList{10}.piston = 150e-9
% DM.segList{20}.tiptilt = [1 2]/2000/NN;

JWST.touch;
F.planewave*JWST*A;
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
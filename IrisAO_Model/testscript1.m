%% Clear Workspace
clear all;
clc;
close all;

%% Set Initial Parameters and Flags
segpitch = 606e-6;
magnification = 1;

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
mask = DM.segList{1}.Segment.grid_;

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
% DM = IrisAOPTT(DM,1:37,zeros(37,1),zeros(37,1),zeros(37,1),false);
% DM.show;
% DM = IrisAOPTT(DM,1:37,zeros(37,1),ones(37,1)*-10^-3,ones(37,1)*10^-3,true);
% figure; DM.show;
DM = IrisAOPTT(DM,1:37,randn(37,1)*10^-3,randn(37,1)*10^-3,randn(37,1)*10^-3,false);
% figure; DM.show;
% DM.grid;

%% Try to do tip tilt separately from AOSegment (works for magnification=1 only!)
segwidthx = 5.757e-4;
segwidthy = 6.85e-4;
for ii = 1:length(DM.segList)
    % Get angles and pistons
    tip = DM.segList{ii}.tiptilt(1);
    tilt = DM.segList{ii}.tiptilt(2);
    piston = DM.segList{ii}.piston;
    if piston == 0
        piston = 1e-10;
    end
    if tip == 0
        tip = 1e-10;
    end
    if tilt == 0
        tilt = 1e-10;
    end
    % compute heights
    tipheightmax = (segwidthy * sin(tip))/2;
    tiltheightmax = (segwidthx * sin(tilt))/2;
    % make vectors
    tipvec = linspace(-tipheightmax,tipheightmax,floor(segwidthy/DM.dx)+2);
    tiltvec = linspace(-tiltheightmax,tiltheightmax,floor(segwidthx/DM.dx)+12);
    %add zeros to vectors
    sizing = DM.segList{ii}.Segment.size;
    tipdif = (sizing(1))-length(tipvec);
    tiltdif = (sizing(2))-length(tiltvec);
    zerotip1 = zeros(1,tipdif/2);
    zerotilt1 = zeros(1,tiltdif/2);
    tipvec = horzcat(horzcat(zerotip1,tipvec),zerotip1);
    tiltvec = horzcat(horzcat(zerotilt1,tiltvec),zerotilt1);   
    [X,TIP] = meshgrid(tipvec);
    [TILT,YY] = meshgrid(tiltvec);
    % Add in Piston
    OPL = 2*piston;
    nugrid = (TIP+TILT+OPL).*mask;
    grid(DM.segList{ii}.Segment,nugrid);
    
    figure(1)
%     subplot(1,3,1)
%     imagesc(TIP .* DM.segList{ii}.Segment.grid_)
%     sqar;
%     subplot(1,3,2)
%     imagesc(TILT.*DM.segList{ii}.Segment.grid_);
%     sqar;
%     title(sprintf('Segment # %d\n',ii));
%     subplot(1,3,3)
    imagesc(DM.segList{ii}.Segment.grid_);
    title(sprintf('Segment # %d\n',ii));
    sqar;
    axis xy;
    drawnow;
    input('Continue');
end

%% Make an AOGrid to show how it renders
M = AOGrid(DM);
M.spacing(DM.spacing);
M.zero;
for ii = 1:length(DM.segList)
    N = AOGrid(1);
    DM.segList{ii}.Segment.Offset = DM.segList{ii}.Offset;
    N.spacing(DM.segList{ii}.Segment.spacing);
    N.Offset = DM.segList{ii}.Segment.Offset;
    N.grid(DM.segList{ii}.Segment.grid_);
%     N.show;
%     drawnow;
    
    M = M + N;
    M.show;
    drawnow;
    sqar;
    input('Continue');
end
DM.combined.grid(M.grid);
F.planewave * DM * A;
figure;
F.show
colormap(gray);

else
    F_Scal * DM * A;
    F_Scal.show;
end

%% Make a PSF
% THld = lambda/D * 206265; % Lambda/D in arcsecs.
% FOV = 50*THld; % FOV for PSF computation
% PLATE_SCALE = THld/5; % Pixel Size for PSF computation
% [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
% PSFmax = max(PSF(:));
% 
% if verbose == true
%     figure(2)
%     colormap(gray);
% end
% imagesc(thx,thy,log10(PSF/PSFmax),[-4,0]);
% % imagesc(thx,thy,PSF/PSFmax);
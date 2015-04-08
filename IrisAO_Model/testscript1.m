%% Clear Workspace
clear all;
clc;
close all;

%% Set Initial Parameters and Flags
segpitch = 606e-6;
magnification = 1;

verbose = false;
Scalloped_field = false;
DEBUG = false;
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
    % Get something to put on mirror
load System_Tip_Correction_PTTpos.mat;
PTTpos = PTTPositionArray;
% DM = IrisAOPTT(DM,1:37,PTTpos(:,1)*10^-6,PTTpos(:,2)*10^-3,PTTpos(:,3)*10^-3,false);
% DM = IrisAOPTT(DM,1:37,zeros(37,1)*10^-6,zeros(37,1)*-10^-3,zeros(37,1)*10^-3,false);
DM = IrisAOPTT(DM,1:37,zeros(37,1)*10^-6,ones(37,1)*-10^-3,ones(37,1)*-10^-3,false);
% DM = IrisAOPTT(DM,1:37,randn(37,1)*10^-6,randn(37,1)*-10^-3,randn(37,1)*10^-3,true);

%% Do piston tip tilt separately from AOSegment ()
%Segment widths in the model (for reference, segwidthy should be 7.00e-4 meters)
segwidthx = 5.757e-4 * magnification;
segwidthy = 6.85e-4 * magnification;
%Segment Coordinate System
[xseg,yseg] = DM.segList{1}.Segment.coords;
%Build an AOGrid Object to do the rendering
M = AOGrid(DM);
M.name = 'IrisAO Mirror PTT';
M.spacing(DM.spacing);
M.zero;
%Go segment by segment.....
for ii = 1:length(DM.segList)
    % Get angles and pistons
    tip = DM.segList{ii}.tiptilt(1);
    tilt = DM.segList{ii}.tiptilt(2);
    piston = DM.segList{ii}.piston;
    
% if mirror is set to flat, add a little piston so it keeps the 
% segments (sort of a bug, but it has to be done until a better fix is found)
    if piston == 0
        if tip == 0
            if tilt == 0
                piston = 1e-10; %1 Angstrom of piston
            end
        end
    end

    % compute edge displacements
    tipheightmax = (segwidthy * sin(tip))/2;
    tiltheightmax = (segwidthx * sin(tilt))/2;
    
    % make vectors for displacement, across only defined segment positions
    %(number of points is tweaked manually to ensure all of the smoothed edge segments are included
    tipvec = linspace(-tipheightmax,tipheightmax,floor(segwidthy/DM.dx)+2);
    tiltvec = linspace(-tiltheightmax,tiltheightmax,floor(segwidthx/DM.dx)+12);
    
    %add zeros to vectors to bring them up to the length of the matrix sides
    sizing = DM.segList{ii}.Segment.size;
    tipdif = (sizing(1))-length(tipvec);
    tiltdif = (sizing(2))-length(tiltvec);
    zerotip1 = zeros(1,tipdif/2);
    zerotilt1 = zeros(1,tiltdif/2);
    tipvec = horzcat(horzcat(zerotip1,tipvec),zerotip1);
    tiltvec = horzcat(horzcat(zerotilt1,tiltvec),zerotilt1); 
    
    % Compute Tip/Tilt Grids
    [X,TIP] = meshgrid(tipvec);
    [TILT,Y] = meshgrid(tiltvec);
    
    % Compute Piston
    OPL = 2*piston; %mirror OPL
    k = (2*pi)/lambda;
    % Calculate Correct Segment Shape and Apply
    nugrid = exp((1i*k).*(TIP+TILT+OPL)).*mask;
    grid(DM.segList{ii}.Segment,nugrid);
    
    % The nasty bit of rendering the mirror
    N = AOGrid(1);
    % Grab the offsets so the code knows where to put the segment
    DM.segList{ii}.Segment.Offset = DM.segList{ii}.Offset;
    % Ensure the spacing is the same
    N.spacing(DM.segList{ii}.Segment.spacing);
    % Apply the offset
    N.Offset = DM.segList{ii}.Segment.Offset;
    % Set the Segment Shape
    N.grid(DM.segList{ii}.Segment.grid_);
    % Add the Segment to the Mirror
    M = M + N;
    
    % A snippet that allows the visual checking of what exactly is being
    % set to all the .grid_ properties.  This is where a lot of the error
    % was happening in the previous model, other than the computation that
    % does the projected piston thing John talks about. If the first 3
    % plots are the same, and the most recent segment added to the 4th
    % frame is the same as those 3, the rendering is working correctly
if DEBUG == true
    figure(1)
    subplot(2,2,1)
    imagesc(xseg,yseg,abs(nugrid));
    sqar;
    axis xy;
    colorbar;
    caxis([-4e-6,4e-6]);
    title(sprintf('nugrid for \nSegment # %d\n',ii));
    subplot(2,2,2)
    imagesc(xseg,yseg,abs(DM.segList{ii}.Segment.grid_));
    sqar;
    axis xy;
    colorbar;
    caxis([-4e-6,4e-6]);
    title(sprintf('Segment.grid_ for \nSegment # %d\n',ii));
    subplot(2,2,3);
    imagesc(xseg,yseg,abs(N.grid));
    sqar;
    axis xy;
    colorbar;
    caxis([-4e-6,4e-6]);
    title(sprintf('N.grid_ for \nSegment # %d\n',ii));
    subplot(2,2,4)
    imagesc(x,y,abs(M.grid));
    sqar;
    axis xy;
    title(sprintf('Adding \nSegment # %d\n',ii));
    colorbar;
    caxis([-4e-6,4e-6]);
    drawnow;
    pause(0.5);
else
    % The rendering is being trusted by the user. Just plot the Segment
    % shape and the segments as they are added to the mirror
    if verbose == true
        figure(1)
        subplot(1,2,1)
        imagesc(xseg,yseg,abs(DM.segList{ii}.Segment.grid_));
        sqar;
        axis xy;
        title(sprintf('PTT for \nSegment # %d\n',ii));
        subplot(1,2,2)
        imagesc(x,y,abs(M.grid));
        sqar;
        axis xy;
        title(sprintf('Adding Segment #%d\nto IrisAO Mirror\n',ii));
        drawnow;
        pause(0.1);
    end
end
    
end

% Apply the rendered mirror shape to the DM variable
DM.combined.grid(M.grid);

F.planewave * DM * A;
if verbose == true
    figure(2);
else
    figure(1);
end

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
    figure(3)
    colormap(gray);
else
    figure(2)
end

% imagesc(thx,thy,PSF/PSFmax);
imagesc(thx,thy,log10(PSF/PSFmax),[-4,0]);
colormap(gray);
axis xy;
sqar;

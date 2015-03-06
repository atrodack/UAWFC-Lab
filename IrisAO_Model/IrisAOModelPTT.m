function Positioned_Segments = IrisAOModelPTT(Segment,SegPitch,Piston,Tip,Tilt,gain)
% Positioned_Segments = IrisAOModelPTT(Segment,Piston,Tip,Tilt)
% 
% Author: Alexander Rodack
% 
% Date of Creation: 10/3/2014
% 
% Description:
% This function takes values for tip, tilt in micron, mrad, mrad,
% and applies them to the IrisAO DM model segment masks.
% 
% Input Parameters
% Segment: This is the Segment Mask to apply tip/tilt to created by the 
% IrisAO Model.
% 
% SegPitch: This is the Diameter of the Mirror Segment
%
% Piston: Length in micron to apply for Piston
% 
% Tip: Angle in mrad to apply for Tip
% 
% Tilt: Angle in mrad to apply for Tilt
% 
% gain: Set gain value to apply to make motions more visible. NOTE: only 
% change this from 1 when using this code to visualize motions.  If used in
% simulation, leave at 1 to maintain correct calculated values.
% 
% Function Output
% Positioned_Segments: A cell array containing the positioned segment masks
% 
% Version: 0.1.5




% *************************************************************************

if nargin < 6
    gain = 1;
end

Piston = Piston * 10^-6;
Tip = Tip *10^-3;
Tilt = Tilt * 10^-3;

% *************************************************************************
% Set Segment Radii
yradius = SegPitch/2; %m
xradius = SegPitch/2; %m

% *************************************************************************
% Obtain Parameters
numrows = size(Segment,1);
numcols = size(Segment,2);
testmaskpiston = Segment;

% This is determining how many rows/cols there are 1's in.  This will be
% utilized in the new algorithm of applying the correct tip/tilt to the
% segments
testmaskrows = Segment;
counterrows = 0;
tol = 0;
tol2 = 0;

for k = 1:numrows
    testvectorrows = testmaskrows(k,:);
    if max(testvectorrows) > tol2
        counterrows = counterrows+1;
    end
end
% fprintf('The number of rows with Segment Present is %d \n',counterrows);

testmaskcols = Segment;
countercols = 0;
for k = 1:numcols
    testvectorcols = testmaskcols(:,k);
    if max(testvectorcols) > tol2
        countercols = countercols+1;
    end
end
% fprintf('The number of cols with Segment Present is %d \n',countercols);

% *************************************************************************
% Calculate Pistons Needed for Tip
Piston_Tip = yradius * tan(Tip);
% fprintf('Segment edge must be pistoned by %0.5f micron to achieve Tip of %0.2f mrad\n',Piston_Tip,Tip);
tipvec = linspace(-Piston_Tip,Piston_Tip,countercols);
% tipvec = tipvec + 0.5;
% *************************************************************************
% Calculate Pistons Needed for Tilt
Piston_Tilt = xradius * tan(Tilt);
% fprintf('Segment edge must be pistoned by %0.5f micron to achieve Tilt of %0.2f mrad\n',Piston_Tilt,Tilt);
tiltvec = linspace(-Piston_Tilt,Piston_Tilt,counterrows);
% tiltvec = tiltvec + 0.5;

% *************************************************************************
% Apply to Segment
% Tilt
index = 1;
for k = 1:numrows
    testvectorrows = testmaskrows(k,:);
    if max(testvectorrows) > tol
        for jj = 1:numcols
            if testmaskrows(k,jj) > tol2
                testmaskrows(k,jj) = tiltvec(index);
            elseif testmaskrows(k,jj) > tol2
%                 testmaskrows(k,jj) = tiltvec(index)+0.5*testmaskrows(k,jj);
            end
        end
        index = index+1;
    end
end
% Tip
index = 1;
for k = 1:numcols
    testvectorcols = testmaskcols(:,k);
    if max(testvectorcols) > tol
        for jj = 1:numrows
            if testmaskcols(jj,k) > tol2
                testmaskcols(jj,k) = tipvec(index);
            elseif testmaskcols(jj,k) > tol2
%                 testmaskcols(k,jj) = tipvec(index) + 0.5*testmaskcols(k,jj);
            end
        end
        index = index+1;
    end
end

% Combine into one Tip/Tilt mask
TTmask = testmaskcols + testmaskrows;

% Piston
index = 1;
for k = 1:numrows
    testvectorrows = testmaskpiston(k,:);
    if max(testvectorrows) > tol
        for jj = 1:numcols
            if testmaskpiston(k,jj) > tol2
                testmaskpiston(k,jj) = Piston;
            end
        end
        index = index+1;
    end
end
Pmask = testmaskpiston;

% Create Final Segment Position Mask with Piston/Tip/Tilt
PTTmask = TTmask + Pmask;
% fprintf('Segment Complete\n\n');

% ************************************************************************* 
% Create Function Output and Apply Gain
Positioned_Segments = gain*PTTmask; %apply gain to see movement


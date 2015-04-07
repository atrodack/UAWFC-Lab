function [ DM ] = IrisAOPTT( DM, segList, pistonList, tipList, tiltList, bump)
%[ DM ] = IrisAOPTT( DM, segList, pistonList, tipList, tiltList )
%   Applys Piston, Tip, and Tilt to an IrisAO model. This will do it in
%   hardware order, not in the Model's order. Input segList can be 1:37,
%   and this will map to the corresponding segment in DM.segList.
%
% pistonList is in units of meters---> on order of microns
% tipList and tiltList are in units of radians ----> on order of mrad

%% Check Inputs
if nargin == 1
    segList = 1:37;
    pistonList = zeros(37,1);
    tipList = zeros(37,1);
    tiltList = zeros(37,1);
    bump = false;
    fprintf('Flattening the Mirror\n');
elseif nargin == 3
    tipList = zeros(37,1);
    tiltList = zeros(37,1);
    bump = false;
    fprintf('Applying Piston Only\n');
elseif nargin == 4
    tiltList = zeros(37,1);
    bump = false;
    fprintf('Applying Piston and Tip\n');
elseif nargin == 5
    fprintf('Applying Piston, Tip, and Tilt\n');
    bump = false;
elseif nargin == 6
    if bump == true
        fprintf('Bumping Segments\n');
    else
        fprintf('Applying Piston, Tip, and Tilt\n');
    end
else
    error('Number of input arguments is incorrect');
end

if ~isa(DM,'AOAperture')
    error('DM must be AOSim2 class AOAperture');
end

%% Load in the Mapping Vector
% load('IrisAO_MAP.mat');
load('IrisAO_STUFF.mat');
figure;
%% Do the Mapping and Apply PTT
if bump == false
    jj = 1;
    for ii = 1:37
        while(jj <= 37)
            if IrisAO_MAP(jj) == segList(ii)
                mapped_segment = jj;
                jj = 50;
            else
                jj = jj + 1;
            end
        end
        
        DM.segList{ii}.piston = pistonList(mapped_segment);
        DM.segList{ii}.tiptilt = [tipList(mapped_segment),tiltList(mapped_segment)];
        jj = 1;
    end
else
    jj = 1;
    for ii = 1:37
        while(jj <= 37)
            if IrisAO_MAP(jj) == segList(ii)
                mapped_segment = jj;
                jj = 50;
            else
                jj = jj + 1;
            end
        end
        
        DM.segList{ii}.piston = DM.segList{ii}.piston + pistonList(mapped_segment);
        DM.segList{ii}.tiptilt = DM.segList{ii}.tiptilt + [tipList(mapped_segment),tiltList(mapped_segment)];
        jj = 1;
    end
end

%% Go Ahead and Render It
DM.touch;
DM.render;

end


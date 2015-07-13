function [PTT_mapped] = mapSegments(PTTpos, numRings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[sy,sx] = size(PTTpos);
PTT_mapped = zeros(sy,sx);
if nargin == 1
%     numRings = 3; %default to 37 segments
    numsegs = length(PTTpos);
    switch numsegs
        case 1
            numRings = 0;
            
        case 7
            numRings = 1;
            
        case 19
            numRings = 2;
            
        case 37
            numRings = 3;
    end
end
% Load in Mapping Data
%         load('IrisAO_SegMap.mat');
switch numRings
    case 0
        IrisAO_SegMap = [1];
        
    case 1
        IrisAO_SegMap = [6 7 5 1 2 4 3]; % 7 segment
        
    case 2
        IrisAO_SegMap = [16 17 18 15 6 7 19 14 5 1 2 8 13 4 3 9 12 11 10]; % 19 segment
        
    case 3
        IrisAO_SegMap = [32 33 34 35 31 16 17 18 36 30 15 6 7 19 37 29 14 5 1 2 8 20 28 13 4 3 9 21 27 12 11 10 22 26 25 24 23]; %37 segment
        
    otherwise
        error('This is not the mirror you are looking for...');
        
end
% Map the PTT matrix from hardware to software order
numSeg = length(IrisAO_SegMap);
for ii = 1:numSeg
    mapped_segment = IrisAO_SegMap(ii);
    PTT_mapped(ii,1:3) = PTTpos(mapped_segment,:);
end

end


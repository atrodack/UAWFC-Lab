function [ PTT_mapped ] = mapSegments( PTTpos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[sy,sx] = size(PTTpos);
PTT_mapped = zeros(sy,sx);

 % Load in Mapping Data
%         load('IrisAO_SegMap.mat');
        IrisAO_SegMap = [32 33 34 35 31 16 17 18 36 30 15 6 7 19 37 29 14 5 1 2 8 20 28 13 4 3 9 21 27 12 11 10 22 26 25 24 23];
        % Map the PTT matrix from hardware to software order
        for ii = 1:37
            mapped_segment = IrisAO_SegMap(ii);
            PTT_mapped(ii,1:3) = PTTpos(mapped_segment,:);
        end

end

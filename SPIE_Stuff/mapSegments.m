function [ PTT_mapped ] = mapSegments( PTTpos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[sy,sx] = size(PTTpos);
PTT_mapped = zeros(sy,sx);

 % Load in Mapping Data
        load('IrisAO_SegMap.mat');
        % Map the PTT matrix from hardware to software order
        for ii = 1:37
            mapped_segment = IrisAO_SegMap(ii);
            PTT_mapped(ii,1:3) = PTTpos(mapped_segment,:);
        end

end


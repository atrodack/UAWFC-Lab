function [ AOdOTF ] = loadinstored( AOdOTF,cell_num )
%LOADINSTORED Summary of this function goes here
%   Detailed explanation goes here


AOdOTF.PSF0 = AOdOTF.storedPSF0{cell_num};
AOdOTF.PSF1 = AOdOTF.storedPSF1{cell_num};
AOdOTF.usedata = true;
AOdOTF.useData;


end


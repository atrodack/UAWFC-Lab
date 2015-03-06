function [ PTTpos ] = IrisAOComputeZernPositions( wavelength, Zernike_Nolls, Zernike_Coefficient_waves )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lambda = wavelength;
ZernModes = Zernike_Nolls;
ZernCoeffs = Zernike_Coefficient_waves;

NumSegments = 37;

    % Initalize PTTpos
PTTpos = zeros(NumSegments,3);

    % Initialize Coefficient Array (needs to be 21x1)
    ZernikeCoefficientArray = zeros(21,1);
    
    % Convert Coefficients from waves to distances
    ZernCoeffs = ZernCoeffs*lambda;
    
    % Add in Coefficients
    ZernikeCoefficientArray(ZernModes) = ZernCoeffs;

    % Load in Zernike Basis Data and Get PTT Positions
    PTTPositionArray = GetPTT489ZernikePostions(ZernikeCoefficientArray);
    
    % Shrink to the size of PTT111
    PTTPositionArray = PTTPositionArray(1:37,:);
    PTTpos = PTTPositionArray;

    fprintf('*****************************************************************\n');
    display('Segment Positions Calculated');
    fprintf('*****************************************************************\n');

end


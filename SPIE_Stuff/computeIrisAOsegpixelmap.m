function [ pixel_seg_map, Areal_Averaging_radius ] = computeIrisAOsegpixelmap( DM, Aperture, pokeseg,numsegs)
% [ pixel_seg_map,Areal_Averaging_radius ] = computeIrisAOsegpixelmap( DM, Aperture, pokeseg)
%   function for computing the Segment-pixel grid for AOSim2 IrisAO Model
%   in dOTF

if nargin < 4
    numsegs = 37;
end

SPACING = DM.spacing;
lambda = AOField.HeNe_Laser;


F1 = AOField(2^11);
F1.spacing(SPACING);
F1.FFTSize = 2^11;
F1.lambda = lambda;

F2 = F1.copy;


% pixel_seg_map = cell(length(DM.segList),2);
pixel_seg_map = cell(length(DM.segList),1);

for n = 1:numsegs
    PTTpos_flat = zeros(length(DM.segList),3);
    PTTpos_poked = zeros(length(DM.segList),3);
%     PTTpos_poked(pokeseg,1) = lambda/4;
PTTpos_poked(pokeseg,3) = 1e-3;
    fprintf('Segment # %d Poked\n',n);
    
    if n ~= pokeseg
        PTTpos_flat(n,1) = 0.5;
        PTTpos_poked(n,1) = 0.5;
        PTT_flat = mapSegments(PTTpos_flat);
        PTT_poked = mapSegments(PTTpos_poked);
        
        DM.setIrisAO(PTT_flat);
        
        F1.planewave * Aperture * DM;
        grid1 = F1.grid;
        PSF1 = abs(fftshift(fft2(fftshift(grid1)))).^2;
        
        DM.setIrisAO(PTT_poked);
        
        F2.planewave * Aperture * DM;
        grid2 = F2.grid;
        PSF2 = abs(fftshift(fft2(fftshift(grid2)))).^2;
        
        OTF1 = fftshift(fft2(fftshift(PSF1)));
        OTF2 = fftshift(fft2(fftshift(PSF2)));
        
        dOTF = OTF1 - OTF2;
%         dOTF = -1i * conj(dOTF);
        dOTF = conj(dOTF);
        
        figure(1);
        plotComplex(dOTF,6);
        axis xy;
        
        if n > 1
        hold on;
        for m = 1:n-1
            if m ~= pokeseg
                plot(pixel_seg_map{m}(2),pixel_seg_map{m}(1),'r*');
            end
        end
        hold off;
        end
        
        drawnow;
        
        display('Pick the center point of the Segment in the leftmost pupil');
        point_left = pickPoint;
%         display('Pick the center point of the Segment in the rightmost pupil');
%         point_right = pickPoint;
        
        if n == 1
            display('Pick the point at the edge of the Segment');
            edge_point = pickPoint;
            Areal_Averaging_radius = sqrt((edge_point(2) - point_left(2))^2 + (edge_point(1) - point_left(1))^2);
        end
        
        pixel_seg_map{n,1} = point_left;
%         pixel_seg_map{n,2} = point_right;
        
    else
        pixel_seg_map{n,1} = [];
%         pixel_seg_map{n,2} = [];
        
        
        
    end
    
    
end


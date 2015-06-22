function [ pixel_seg_map, Areal_Averaging_radius ] = computeIrisAOsegpixelmap( DM, Aperture, pokeseg, FOV, PLATE_SCALE, FoV)
% [ pixel_seg_map ] = computeIrisAOsegpixelmap( DM, pokeseg)
%   function for computing the Segment-pixel grid for AOSim2 IrisAO Model
%   in dOTF
SPACING = DM.spacing;
lambda = AOField.HeNe_Laser;


F = AOField(2^12);
F.spacing(SPACING);
F.FFTSize = 2^12;
F.lambda = lambda;

F2 = F.copy;


% pixel_seg_map = cell(length(DM.segList),2);
pixel_seg_map = cell(length(DM.segList),1);

for n = 1:length(DM.segList)
    PTTpos_flat = zeros(length(DM.segList),3);
    PTTpos_poked = zeros(length(DM.segList),3);
    PTTpos_poked(pokeseg,1) = 0.5;
    fprintf('Segment # %d Poked\n',n);
    
    if n ~= pokeseg
        PTTpos_flat(n,1) = 0.5;
        PTTpos_poked(n,1) = 0.5;
        PTT_flat = mapSegments(PTTpos_flat);
        PTT_poked = mapSegments(PTTpos_poked);
        
        DM.PTT(PTT_flat);
        DM.touch;
        DM.render;
        
        F.planewave * Aperture * DM;
        PSF1 = F.mkPSF(FOV,PLATE_SCALE);
        PSF1max = max(max(PSF1));
        % PSF1 = PSF1 / PSF1max;
        PSF1plot = log10(PSF1/PSF1max);
        
        DM.PTT(PTT_poked);
        DM.touch;
        DM.render;
        
        F2.planewave * Aperture * DM;
        PSF2 = F2.mkPSF(FOV,PLATE_SCALE);
        PSF2max = max(max(PSF2));
        % PSF2 = PSF2 / PSF2max;
        PSF2plot = log10(PSF2/PSF2max);
        
        F.touch; F2.touch;
        F.grid(PSF1); F2.grid(PSF2);
        
        OTF1 = F.mkOTF2(FoV,SPACING(1));
        OTF2 = F2.mkOTF2(FoV,SPACING(1));
        
        F.touch; F2.touch;
        
        dOTF = OTF1 - OTF2;
        
        figure(1);
        plotComplex(dOTF,3);
        
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

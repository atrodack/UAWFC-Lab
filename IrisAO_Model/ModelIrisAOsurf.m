function [  ] = ModelIrisAOsurf( DM,PTTpos )
%[  ] = ModelIrisAOsurf( PTTpos )
%   This function will take the PTT positions of the DM segments and apply
%   them to the mirror model. The output of this function is a subplot of
%   the surface profile of the DM from 2 perspectives


if ~isa(DM,'AOAperture')
    error('DM needs to be an IrisAO Model');
end

%% Apply Tip/Tilt to Individual Segment Masks
% Set Visualization Gain
gain = 1;
% Initialize cell array to hold positioned segments
Positioned_Segments = cell(1,37);

% if ~abs(PTTpos(SegNum,1)) > 0
%     PTTpos(SegNum,1) = 1;
% end
IrisAONumbering; %gets IrisAO_MAP
for k = 1:37
% Segment
%     SegNum = IrisAO_MAP(k);
    Segment = DM.segList{k}.Segment.grid;
    
% Tip angle in mrad
    Tip = PTTpos(SegNum,2);
    
% Tilt angle in mrad
    Tilt = PTTpos(SegNum,3);
    
% Piston in micron
    Piston = PTTpos(SegNum,1);
    
% Apply Piston/Tip/Tilt
    Positioned_Segment = IrisAOModelTT(Segment,Piston,Tip,Tilt,gain);
    Positioned_Segments{k} = Positioned_Segment;
end

%% Create Final Mask
% This mask contains the information of all 37 segments after they have
% been positioned.

FinalMask = Positioned_Segments{1};
for ii = 2:37
    FinalMask = FinalMask + Positioned_Segments{ii};
end

%% Display Mirror Surface
% figure(10)
% % Plot Profile View
% subplot(1,2,1)
imagesc(x,x,FinalMask);
axis xy
hold on
plot(SegBoundaries(:,1),SegBoundaries(:,2),'ks');
hold off;
daspect([1 1 1]);
title(sprintf('Profile View with Z values scaled by %0.2f',gain));
% 
% % Plot DM Surface
% subplot(1,2,2)
% surf(FinalMask,'LineStyle','none');
% zlim([-6 6]);
% daspect([1 1 1]);
% lt=light();
% set(lt,'Position',[0,0,10]);
% title(sprintf('DM Shape with Z values scaled by %0.2f \n\n',gain));
% % colormap(gray);



% mirror_shape = FinalMask;
% phi = mirror_shape;
% A = 1;
% field = A .* exp(2*pi*1i.*phi);
% 
% field = field .* mirror_shape;
% field = padarray(field,[1024,1024]);
% PSF = abs(fftshift(fft2(fftshift(field)))).^2;
% maxval = max(PSF(:));
% PSF = PSF./maxval;
% figure;
% imagesc(PSF.^0.5);
% axis off;
% axis xy;
% colormap(jet)

end


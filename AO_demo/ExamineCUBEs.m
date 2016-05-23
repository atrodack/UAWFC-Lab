function [CUBE1, CUBE2, dOTF] = ExamineCUBEs(CUBE1, CUBE2, section_num)

nframes = size(CUBE1.PSFs,4);
dOTF = zeros(256,256);
data = 0;
% OTF_CUBE1 = zeros(256,256,1,nframes);
% OTF_CUBE2 = zeros(256,256,1);

if section_num == 1
    %% display PSFs
    figure(20);
    for n = 1:nframes
        imagesc(CUBE1.PSFs(:,:,1,n));
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
    
    for n = 1:nframes
        imagesc(CUBE2.PSFs(:,:,1,n));
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
elseif section_num == 2
    %% Compute OTFs
    for n = 1:nframes
        
        PSF1 = double(CUBE1.PSFs(:,:,1,n));
        PSF2 = double(CUBE2.PSFs(:,:,1,n));
        
        % HeNe (needs recalibration after movement of flip mirror/OD filter)
        % centerpoint = [272,305]; %[y,x] in ij
        
        % SuperK (needs recalibration after input of OD filter)
        centerpoint = [254,316]; %no OD filter
        
        % Cenering PSFs and cropping to 256x256
        PSF1 = PSF1(centerpoint(1)-127:centerpoint(1)+128,centerpoint(2)-127:centerpoint(2)+128);
        PSF2 = PSF2(centerpoint(1)-127:centerpoint(1)+128,centerpoint(2)-127:centerpoint(2)+128);
        
        CUBE1.cropped_PSFs(:,:,1,n) = PSF1;
        CUBE2.cropped_PSFs(:,:,1,n) = PSF2;
        
        % Shift PSFs to Corners (fftshift by hand to make sure pixel allignment is correct)
        centerpoint2 = [129,129];
        PSF1 = circshift(PSF1,1-centerpoint2);
        PSF2 = circshift(PSF2,1-centerpoint2);
        
        % Compute the OTFs and zero central pixel
        OTF1 = fftshift(fft2(PSF1));
        OTF2 = fftshift(fft2(PSF2));
        
        OTF1(centerpoint2(1),centerpoint2(2)) = 0;
        OTF2(centerpoint2(1),centerpoint2(2)) = 0;
        
        CUBE1.OTFs(:,:,1,n) = OTF1;
        CUBE2.OTFs(:,:,1,n) = OTF2;
    end
    
elseif section_num == 3
    %% Show OTFs
    figure(20)
    for n = 1:nframes
        plotComplex(CUBE1.OTFs(:,:,1,n),4);
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
    
    
    
    for n = 1:nframes
        plotComplex(CUBE2.OTFs(:,:,1,n),4);
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
    
elseif section_num == 4
    %% Average PSF and OTF
    PSF1_avg = 0;
    PSF2_avg = 0;
    OTF1_avg = 0;
    OTF2_avg = 0;
    
    for n = 1:nframes
        PSF1_avg = PSF1_avg + double(CUBE1.cropped_PSFs(:,:,1,n));
        PSF2_avg = PSF2_avg + double(CUBE2.cropped_PSFs(:,:,1,n));
        OTF1_avg = OTF1_avg + CUBE1.OTFs(:,:,1,n);
        OTF2_avg = OTF2_avg + CUBE2.OTFs(:,:,1,n);
    end
    
    PSF1_avg = PSF1_avg / nframes;
    OTF1_avg = OTF1_avg / nframes;
    CUBE1.PSF = PSF1_avg;
    CUBE1.OTF = OTF1_avg;
    
    PSF2_avg = PSF2_avg / nframes;
    OTF2_avg = OTF2_avg / nframes;
    CUBE2.PSF = PSF2_avg;
    CUBE2.OTF = OTF2_avg;
    
elseif section_num == 5
    %% Compute dOTF
    % maxval1 = mean(mean(real(CUBE1.OTF)));
    % maxval2 = mean(mean(real(CUBE2.OTF)));
    %
    % if maxval1 > maxval2
    %     scale = maxval1 / maxval2
    %     CUBE2.OTF = scale * CUBE2.OTF;
    % else
    %     scale = maxval2 / maxval1
    %     CUBE1.OTF = scale * CUBE1.OTF;
    % end
    
    dOTF = CUBE2.OTF - CUBE1.OTF;
    dOTF(129,:) = 0;
    dOTF(:,129) = 0;
    dOTF = -1i * conj(dOTF);
    figure(21);
    plotComplex(dOTF,3)
    
    
elseif section_num == 6
    %% Shift OTF Stacks
    CEN = [129,129];
    figure(20)
    plotComplex(CUBE1.OTFs(:,:,1,1),6);
    axis xy;
    title(sprintf('Pick Point\nat radius of OTF'));
    PT = pickPoint(1);
    
    RADIUS = sqrt((abs(CEN(2) - PT(2))).^2 + (abs(CEN(1) - PT(1))).^2);
    TT1 = zeros(1,2,nframes);
    TT2 = zeros(1,2,nframes);
    for n = 1:nframes
        [CUBE1.OTFs_(:,:,1,n),TT1(:,:,n)] = vernierAlignOTF(CUBE1.OTFs(:,:,1,n),RADIUS,CEN);
        [CUBE2.OTFs_(:,:,1,n),TT2(:,:,n)] = vernierAlignOTF(CUBE2.OTFs(:,:,1,n),RADIUS,CEN);
    end
    
    CUBE1.TT = TT1;
    CUBE2.TT = TT2;
    
elseif section_num == 7
    %% Show Shifted OTFs
    figure(20)
    for n = 1:nframes
        plotComplex(CUBE1.OTFs_(:,:,1,n),4);
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
    
    
    
    for n = 1:nframes
        plotComplex(CUBE2.OTFs_(:,:,1,n),4);
        axis xy;
        sqar;
        title(sprintf('Frame %d',n));
        drawnow;
        pause(0.1);
    end
    
elseif section_num == 8
    %% Average Shifted OTFs
    
    OTF1_avg = 0;
    OTF2_avg = 0;
    
    for n = 1:nframes
        OTF1_avg = OTF1_avg + CUBE1.OTFs_(:,:,1,n);
        OTF2_avg = OTF2_avg + CUBE2.OTFs_(:,:,1,n);
    end
    
    OTF1_avg = OTF1_avg / nframes;
    CUBE1.OTF_ = OTF1_avg;
    
    OTF2_avg = OTF2_avg / nframes;
    CUBE2.OTF_ = OTF2_avg;
    
elseif section_num == 9
    %% Shift Average OTFs
    fprintf('NOT WORKING YET\n');
    % CEN = [129,129];
    % RADIUS = 52;
    % [CUBE1.OTF_,TT1] = vernierAlignOTF(CUBE1.OTF_,RADIUS,CEN);
    % [CUBE2.OTF_,TT2] = vernierAlignOTF(CUBE2.OTF_,RADIUS,CEN);
    
elseif section_num == 10
    %% Compute dOTF
    % maxval1 = mean(mean(real(CUBE1.OTF_)));
    % maxval2 = mean(mean(real(CUBE2.OTF_)));
    
    % if maxval1 > maxval2
    %     scale = maxval1 / maxval2
    %     CUBE2.OTF_ = scale * CUBE2.OTF_;
    % else
    %     scale = maxval2 / maxval1
    %     CUBE1.OTF_ = scale * CUBE1.OTF_;
    % end
    
    dOTF = CUBE2.OTF_ - CUBE1.OTF_;
    dOTF(129,:) = 0;
    dOTF(:,129) = 0;
    dOTF = -1i * conj(dOTF);
    figure(22);
    plotComplex(dOTF,3)
    
elseif section_num == 11
    %% Scan Power
    const = linspace(0.90,1.1,5000);
    x = linspace(1,256,256);
    [X,Y] = meshgrid(x);
    circ = sqrt((X-92).^2 + (Y-126).^2);
    mask = double(circ<10);
    data = zeros(1,length(const));
    GRAB = logical(data);
    for n = 1:length(const)
        patch = (CUBE2.OTF_ - const(n).*CUBE1.OTF_).*mask;
        SELECT = patch(patch~=0);
        data(1,n) = mean(mean(SELECT));
    end
    
    for n = 1:length(data)
        GRAB(n) = isequal(data(n),min(data));
    end
    scale = const(GRAB);
    fprintf('scale = %g\n',scale)
    CUBE1.OTF_ = scale .* CUBE1.OTF_;
    dOTF = CUBE2.OTF_ - CUBE1.OTF_;
    dOTF(129,:) = 0;
    dOTF(:,129) = 0;
    dOTF = -1i * conj(dOTF);
    figure(23);
    plotComplex(dOTF,3)
end %if

end %function
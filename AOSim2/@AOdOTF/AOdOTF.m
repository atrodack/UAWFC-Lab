classdef AOdOTF < AOField
    %AOdOTF A simple addition to AOSim2 that builds a system to use dOTF to
    %compute the wavefront phase.
    
    properties
        finger = [];
        PSF0;
        PSF1;
        OTF0;
        OTF1;
        dOTF;
        Phase;
        OPL;
        A;
        Field;
%         FOV;
%         FoV;
        Plate_Scale;
        thx;
        thy;
        Mask;
        plotMask;
        Mask_interped;
        pupil_center;
        pupil_radius;
        calibration = true;
        usedata = false;
        psf_center;
        storedPSF0 = cell(1,1);
        storedPSF1 = cell(1,1);
        storeddOTF = cell(1,1);
        psfimage_counter = 1;
        dOTFimage_counter = 1;
    end
    
    
    
    methods
        
        % Constructor
        function obj = AOdOTF(A,FOV,Plate_Scale,FoV_OTF,Spacing)
            obj = obj@AOField(A);
            if nargin == 1
                if isa(A,'double')
                    fprintf('You are using captured images for PSFs\n');
                    obj.usedata = true;
%                     obj.psf_center = [260 318];
                elseif isa(A,'AOAperture')
                    obj.A = A;
                    thld = (obj.lambda / (A.segList{1}.Segment.pupils(1,3)) * 206265);
                    obj.FOV =  thld * 100; %make this kind of large to improve dOTF resolution
                    obj.Plate_Scale = thld / 4;
                end
            elseif nargin == 3
                obj.A = A;
                obj.FOV = FOV;
                obj.Plate_Scale = Plate_Scale;
            elseif nargin == 5
                obj.A = A;
                obj.FOV = FOV;
                obj.Plate_Scale = Plate_Scale;
                obj.FoV = FoV_OTF;
                obj.spacing(Spacing);
            else
                error('Incorrect Inputs');
            end
                
        end %Constructor
        
        
        function AOdOTF = setField(AOdOTF,Field)
            AOdOTF.Field = Field;
        end %setField
        
        
        function AOdOTF = calibrateWFS(AOdOTF,y_pos,x_pos,width,Field,ps)
            % This method calbirates the dOTF calculations.  It will create
            % a finger to place in the pupil, and then use a planewave on
            % the pupil to generate the dOTF.  The magnitude of this dOTF
            % is plotted, and used to create the masks needed for phase
            % retrieval. 3 points are selected, the center of the upper
            % pupil, the edge of the upper pupil, and the center of the
            % lower pupil.  The rest is done automatically.  Finally, all
            % of the dOTF calculating properties are cleared.
            
            if nargin == 1
                fprintf('Using Default Values\n');
                [x,y] = AOdOTF.A.coords;
                x_pos = max(x)/2;
                y_pos = max(y)/2;
                width = 0.25;
                Field = AOField(AOdOTF.A);
                ps = 1;
            elseif nargin == 4
                fprintf('Using Default Field\n');
                Field =AOField(AOdOTF.A);
                Field.FFTSize = 1024;
                ps = 1;
            elseif nargin == 5
                ps = 1;
            end
            if AOdOTF.finger == empty
                AOdOTF.create_finger(y_pos,x_pos,width);
            end
            AOdOTF.mkPSF(Field,ps);
            AOdOTF.mkOTF(Field,ps);
            AOdOTF.mkdOTF;
            AOdOTF.plotdOTFframe;
            AOdOTF.mkMask;
            AOdOTF.truephase;
            AOdOTF.calibration = false;
            AOdOTF.cleardOTF;
        end %calibrateWFS
        
        
        function AOdOTF = sense(AOdOTF,Field,ps,ALGO)
            if nargin == 3
                AOdOTF.mkPSF(Field,ps);
                AOdOTF.mkOTF(Field,ps);
                AOdOTF.mkdOTF;
                AOdOTF.truephase;
            elseif nargin == 4
                AOdOTF.mkPSF(Field,ps);
                AOdOTF.mkOTF(Field,ps);
                AOdOTF.mkdOTF;
                AOdOTF.truephase(ALGO);
            end
        end %sense
        
        
        function AOdOTF = calibrateWFS2(AOdOTF,Field)
            fprintf('\n***Calibrating the dOTF WFS***\n');
            F = Field.copy;
            AOdOTF.sense2(F);
            AOdOTF.plotdOTFframe;
            AOdOTF.mkMask;
            AOdOTF.truephase;
            AOdOTF.calibration = false;
            AOdOTF.cleardOTF;
            fprintf('***Calibration Complete***\n');
        end %calibrateWFS2
        
        
        function AOdOTF = precalibratedWFS(AOdOTF,calibration_number)
            fprintf('Using Precalibrated dOTF WFS Object # %d\n',calibration_number);
            switch calibration_number
                case 1
                    load('Calibration_1.mat');
                    AOdOTF.pupil_center = [1552 1005];
                    AOdOTF.pupil_radius = 522.2882;
                    AOdOTF.Mask = MASK;
                    AOdOTF.plotMask = PLOTMASK;
                    AOdOTF.Mask_interped = MASK_INTERPED;
                    AOdOTF.calibration = false;
            end
        end %precalibratedWFS
                    
            
        function AOdOTF = sense2(AOdOTF,Field,Noise,ALGO)
            %AOdOTF = sense2(AOdOTF,Field,Noise,ALGO)
            %Noise is a cell, Noise{1} = flag for using noise, Noise{2} is
            %number of pictures to sum together
            
            if nargin == 2
                SPACING = AOdOTF.spacing;
                AOdOTF.setField(Field);
                Field2 = Field.copy;
                fprintf('Computing PSF\n');
                [AOdOTF.PSF0,AOdOTF.thx,AOdOTF.thy] = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                AOdOTF.setField(Field2 * AOdOTF.finger);
                fprintf('Computing Modified PSF\n');
                AOdOTF.PSF1 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                Field.touch;
                Field.grid(AOdOTF.PSF0);
                Field2.touch;
                Field2.grid(AOdOTF.PSF1);
                AOdOTF.setField(Field);
                fprintf('Computing OTF\n');
                AOdOTF.OTF0 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.setField(Field2);
                fprintf('Computing Modified OTF\n');
                AOdOTF.OTF1 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                fprintf('Computing dOTF\n');
                AOdOTF.mkdOTF;
            elseif nargin == 3
                ALGO = [];
                SPACING = AOdOTF.spacing;
                AOdOTF.setField(Field);
                Field2 = Field.copy;
                fprintf('Computing PSF\n');
                [AOdOTF.PSF0,AOdOTF.thx,AOdOTF.thy] = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                AOdOTF.setField(Field2 * AOdOTF.finger);
                fprintf('Computing Modified PSF\n');
                AOdOTF.PSF1 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                if iscell(Noise)
                    if Noise{1} == true;
                        PSF0_Sum = 0;
                        PSF1_Sum = 0;
                        for n = 1:Noise{2}
                            Noisy_PSF0 = addNoise(AOdOTF.PSF0,Field.grid,true,1,2);
                            Noisy_PSF1 = addNoise(AOdOTF.PSF1,Field.grid,true,1,2);
                            PSF0_Sum = PSF0_Sum + (Noisy_PSF0);
                            PSF1_Sum = PSF1_Sum + (Noisy_PSF1);
                            fprintf('Picture Number %d Complete\n',n);
                        end
                        AOdOTF.PSF0 = abs(PSF0_Sum);
                        AOdOTF.PSF1 = abs(PSF1_Sum);
                    end
                else
                    ALGO = Noise;
                    
                end
                Field.touch;
                Field.grid(AOdOTF.PSF0);
                Field2.touch;
                Field2.grid(AOdOTF.PSF1);
                AOdOTF.setField(Field);
                fprintf('Computing OTF\n');
                AOdOTF.OTF0 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.setField(Field2);
                fprintf('Computing Modified OTF\n');
                AOdOTF.OTF1 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                fprintf('Computing dOTF\n');
                AOdOTF.mkdOTF;
                if ~isempty(ALGO)
                    AOdOTF.truephase(ALGO);
                end
            elseif nargin == 4
                SPACING = AOdOTF.spacing;
                AOdOTF.setField(Field);
                Field2 = Field.copy;
                fprintf('Computing PSF\n');
                AOdOTF.PSF0 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                AOdOTF.setField(Field2 * AOdOTF.finger);
                fprintf('Computing Modified PSF\n');
                AOdOTF.PSF1 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                if Noise{1} == true;
                    PSF0_Sum = 0;
                    PSF1_Sum = 0;
                    for n = 1:Noise{2}
                        Noisy_PSF0 = addNoise(AOdOTF.PSF0,Field.grid,true,1,2.5);
                        Noisy_PSF1 = addNoise(AOdOTF.PSF1,Field.grid,true,1,2.5);
                        PSF0_Sum = PSF0_Sum + (Noisy_PSF0);
                        PSF1_Sum = PSF1_Sum + (Noisy_PSF1);
                        fprintf('Picture Number %d Complete\n',n);
                    end
                    AOdOTF.PSF0 = abs(PSF0_Sum);
                    AOdOTF.PSF1 = abs(PSF1_Sum);
                end
                
                Field.touch;
                Field.grid(AOdOTF.PSF0);
                Field2.touch;
                Field2.grid(AOdOTF.PSF1);
                AOdOTF.setField(Field);
                fprintf('Computing OTF\n');
                AOdOTF.OTF0 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.setField(Field2);
                fprintf('Computing Modified OTF\n');
                AOdOTF.OTF1 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                fprintf('Computing dOTF\n');
                AOdOTF.mkdOTF;
                fprintf('Computing the WFS Phase using the %s Method\n',ALGO);
                AOdOTF.truephase(ALGO);
            end
        end %sense2
        
        
        function AOdOTF = sense_coronagraph(AOdOTF,Field,FPM,Lyot,ALGO)
            if nargin == 4
                SPACING = AOdOTF.spacing;
                Field2 = Field.copy;
                Field = gothroughcoronagraph( Field,FPM,Lyot );
                AOdOTF.setField(Field);
                [AOdOTF.PSF0,AOdOTF.thx,AOdOTF.thy] = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                Field2*AOdOTF.finger;
                Field2 = gothroughcoronagraph( Field2,FPM,Lyot );
                AOdOTF.setField(Field2);
                AOdOTF.PSF1 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                Field.touch;
                Field.grid(AOdOTF.PSF0);
                Field2.touch;
                Field2.grid(AOdOTF.PSF1);
                AOdOTF.setField(Field);
                AOdOTF.OTF0 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.setField(Field2);
                AOdOTF.OTF1 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.mkdOTF;
            elseif nargin == 5
                SPACING = AOdOTF.spacing;
                Field2 = Field.copy;
                Field = gothroughcoronagraph( Field,FPM,Lyot );
                AOdOTF.setField(Field);
                [AOdOTF.PSF0,AOdOTF.thx,AOdOTF.thy] = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                Field2*AOdOTF.finger;
                Field2 = gothroughcoronagraph( Field2,FPM,Lyot );
                AOdOTF.setField(Field2);
                AOdOTF.PSF1 = AOdOTF.Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                Field.touch;
                Field.grid(AOdOTF.PSF0);
                Field2.touch;
                Field2.grid(AOdOTF.PSF1);
                AOdOTF.setField(Field);
                AOdOTF.OTF0 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.setField(Field2);
                AOdOTF.OTF1 = AOdOTF.Field.mkOTF2(AOdOTF.FoV,SPACING(1));
                AOdOTF.mkdOTF;
                AOdOTF.truephase(ALGO);
            end
            
        end %sense_coronagraph
        
        
        function AOdOTF = create_finger(AOdOTF,y_pos,x_pos,width)
            %Creates a pupil finger. y_pos is the height of the finger,
            %x_pos is the x position of the pupil, and width is the width
            %of the finger.
            
            if nargin == 1
                [x,y] = AOdOTF.A.coords;
                x_pos = max(x)/2;
                y_pos = max(y)/2;
                width = 0.25;
            end
            Finger = AOSegment(AOdOTF.A);
            Finger.name = 'Finger';           
            [X,Y] = Finger.COORDS;
            Finger.grid(~(Y<-y_pos & abs(X-x_pos)<width));
            AOdOTF.finger = Finger;
        end %create_finger
        
        
        function AOdOTF = mkPSF(AOdOTF,Field,ps)
            %Makes PSFs and stores them in the properties.  The first uses
            %the unmodified pupil, the second uses the modified pupil.
            %This is only to be called when PSFs are not input by hand (as
            %in actual lab data images). ps must be AOATMO or AOScreen
            %class from AOSim2, and acts as the phase aberration over the
            %pupil. Field must be AOField class from AOSim2.
            
            if nargin == 2
                ps = 1;
            end
%                 if ~isa(ps,'AOAtmo')
%                     if ~isa(ps,'AOScreen')
%                         error('ps must be of type AOAtmo or AOScreen');
%                     end
%                 end
            if ~isa(Field,'AOField')
                error('Field must be of type AOField');
            end
            
            if isempty(AOdOTF.finger)
                AOdOTF.create_finger;
            end
%             Field.FFTSize = 1024;
            Field * ps * AOdOTF.A;
            AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
            Field.touch;
            
            Field * AOdOTF.finger;
            AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
            Field.touch;
        end %mkPSF
                
        
        function AOdOTF = mkOTF(AOdOTF,Field,ps)
            %Makes OTFs and stores them in properties.  If called with no
            %arguments, it simply Fourier Transforms what is stored in the 
            %PSF properties.  If called with one or two arguments, and the
            %PSF properties are empty, it will use Field and/or ps to
            %create and store PSFs, and use those to compute the OTF
            
            if nargin == 1
                ps = 1;
                Field = 1;
            elseif nargin == 2
                if ~isa(Field,'AOField')
                    error('Field must be of type AOField');
                end
                ps = 1;
                if isempty(AOdOTF.finger)
                    AOdOTF.create_finger;
                end
                
                if isempty(AOdOTF.PSF0)
                    Field.FFTSize = 1024;
                    Field.planewave * ps * AOdOTF.A;
                    AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
                
                if isempty(AOdOTF.PSF1)
                    Field * AOdOTF.finger;
                    AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
            elseif nargin == 3
%                 if ~isa(ps,'AOAtmo')  
%                     if ~isa(ps,'AOScreen')
%                         error('ps must be of type AOAtmo or AOScreen');
%                     end
%                 end
                if ~isa(Field,'AOField')
                    error('Field must be of type AOField');
                end
                if isempty(AOdOTF.finger)
                    AOdOTF.create_finger;
                end
                
                if isempty(AOdOTF.PSF0)
                    Field.FFTSize = 1024;
                    Field.planewave * ps * AOdOTF.A;
                    AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
                
                if isempty(AOdOTF.PSF1)
                    Field * AOdOTF.finger;
                    AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
            end
            
            if AOdOTF.usedata == false
                AOdOTF.OTF0 = fftshift(fft2(fftshift(AOdOTF.PSF0)));
                AOdOTF.OTF1 = fftshift(fft2(fftshift(AOdOTF.PSF1)));
            else
                AOdOTF.OTF0 = fftshift(fft2(AOdOTF.PSF0));
                AOdOTF.OTF1 = fftshift(fft2(AOdOTF.PSF1));
%                 AOdOTF.OTF0(240:242,320:322) = 0;
%                 AOdOTF.OTF1(240:242,320:322) = 0;
                AOdOTF.OTF0(241,321) = 0;
                AOdOTF.OTF1(241,321) = 0;
                AOdOTF.decenterpsfs;
            end
        end %mkOTF
        
        
        function AOdOTF = mkdOTF(AOdOTF)
            %Computes the dOTF and stores it as a complex numbered matrix.
            %It will also store the phase of the dOTF in the Phase
            %property.  If running on 64-bit Ubuntu, the unwrapphase method
            %can be uncommented, and the phase stored will be unwrapped
            %using the unwt method (See unwrapphase method)
            
            AOdOTF.dOTF = (AOdOTF.OTF0) - (AOdOTF.OTF1);
            AOdOTF.Phase = angle(AOdOTF.dOTF);
%             AOdOTF.unwrapphase;
        end %mkdOTF
        
        
        function plotPSFframe(AOdOTF)
            imagesc(AOdOTF.PSF0);
            bigtitle('PSF0',12);
            colormap(gray);
            daspect([1,1,1]);
            axis xy;
        end %plotPSFframe
        
        
        function plotdOTFframe(AOdOTF)
            %Plots the magnitude of the dOTF. Most useful in the
            %calibration of the masks, but can be called if desired.
            
            imagesc(abs(AOdOTF.dOTF));
            bigtitle('dOTF Magnitude',12);
            colormap(gray);
            daspect([1 1 1]);
            axis xy;
%             axis off;
        end %plotdOTFframe
        
        
        function AOdOTF = unwrapphase(AOdOTF,ALGO)
            %Unwrapps Phase.  See function uwrap in AOSim2\utils for more
            %information.  Can only be run on Ubuntu 64-bit OS as compiled.
            if nargin == 1
                AOdOTF.Phase = uwrap(AOdOTF.Phase,'gold');
            elseif nargin == 2
                AOdOTF.Phase = uwrap(AOdOTF.Phase,ALGO);
            end
        end %unwrapphase
        
        
        function AOdOTF = truephase(AOdOTF,ALGO)
            %Retrieves the part of the dOTF phase that is important. This
            %uses the mask to grab the upper pupil conjugate, shifts it to
            %the central pixel, clips the matrix to remove surrounding
            %zeros, and then interpolates it up to the size of the DM
            %pupil (assuming AOdOTF.A is the size of the DM).  This might
            %only work for fingers that cause the dOTF to conjugate so
            %there is an upper left pupil and a lower right pupil. Other
            %configurations are UNTESTED.
            if nargin == 1
                ALGO = 'gold';
            end
            
            phase = AOdOTF.Phase;
            mask = AOdOTF.Mask;
            
            shift_point = AOdOTF.pupil_center;
            radius = round(AOdOTF.pupil_radius);
            [sizex,sizey] = size(phase);
            center = [round(sizex/2),round(sizey/2)];
            shift = [shift_point(1) - center(2),shift_point(2) - center(1)];
            phase_ref = phase(center(1),center(2));
            
            if AOdOTF.calibration == false
                fprintf('\nUnwrapping the Phase\n');
                AOdOTF.unwrapphase(ALGO);
                phase = AOdOTF.Phase - phase_ref;
            end
            
            phase = circshift(phase,-shift);
            
            if AOdOTF.calibration == false
                fprintf('Cropping, Masking, and Resizing Computed Phase\n');
            else
                fprintf('Creating Correctly Scaled Mask\n');
            end
            
            % resize to edge of Pupil
            phase = phase(center(1)-radius-0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);
            
            AOdOTF.Phase = phase;
            AOdOTF.resize_phase_to_Pupil;
            
            AOdOTF.Phase = AOdOTF.Phase .* AOdOTF.Mask_interped;
            
            
            if AOdOTF.calibration == false
                lambda = AOdOTF.Field.lambda;
                k = (2*pi) / lambda;
                fprintf('Computing the OPL\n');
                AOdOTF.OPL = AOdOTF.Phase / k;
            end

        end %truephase
        
        
        function AOdOTF = resize_phase_to_Pupil(AOdOTF)
            %Interpolates the stored phase to the size of the DM pupil
            phase = AOdOTF.Phase;
            [sizex,sizey] = size(phase);
            [xx,yy] = AOdOTF.coords;
            x = linspace(min(xx),max(xx),sizex);
            y = x;
            xq = linspace(min(x),max(x),length(xx));
            yq = xq;
            [X,Y] = meshgrid(x,y);
            [Xq,Yq] = meshgrid(xq,yq);
            AOdOTF.Phase = interp2(X,Y,phase,Xq,Yq);
            
        end %resize_phase_to_Pupil

        function AOdOTF = resize_mask_to_Pupil(AOdOTF)
            mask = AOdOTF.Mask;
            
            shift_point = AOdOTF.pupil_center;
            radius = round(AOdOTF.pupil_radius);
            [sizex,sizey] = size(mask);
            center = [round(sizex/2),round(sizey/2)];
            shift = [shift_point(1) - center(2),shift_point(2) - center(1)];
            
            mask = circshift(mask,-shift);
            mask = mask(center(1)-radius - 0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);
            
            AOdOTF.Mask_interped = mask;
            mask = double(AOdOTF.Mask_interped);           
            [sizex,sizey] = size(mask);
            [xx,yy] = AOdOTF.coords;
            x = linspace(min(xx),max(xx),sizex);
            y = x;
            xq = linspace(min(x),max(x),length(xx));
            yq = xq;
            [X,Y] = meshgrid(x,y);
            [Xq,Yq] = meshgrid(xq,yq);
            AOdOTF.Mask_interped = interp2(X,Y,mask,Xq,Yq);
            
        end
        
        
        function points = calibrate_pupil_mask(AOdOTF)
            %A quick method used to grab and store points off of an image.  
            %It uses the points to calculate the radius of the dOTF pupil 
            %and stores this.
            
            fprintf('\n***Please select point at center of the Pupil***\n');
            central_point = pickPoint;
            fprintf('\n***Please select point at edge of the Pupil***\n');
            edge_point = pickPoint;
            fprintf('\n***Please select point at center of the other Pupil***\n\n');
            central_point2 = pickPoint;
            points = cell(1,3);
            points{1} = central_point;
            points{2} = edge_point;
            points{3} = central_point2;
            AOdOTF.pupil_center = central_point;
            AOdOTF.pupil_radius = sqrt((points{2}(1) - points{1}(1))^2 + (points{2}(2) - points{1}(2))^2);
        end %calibrate_pupil_mask
        
        
        function AOdOTF = mkMask(AOdOTF)
            %Creates Boolean logic masks for the dOTF. plotMask is the
            %union of masks for both pupil conjugates to make plotting the
            %unprocessed dOTF phase look nicer.  Mask is the important
            %calculation, and is the piece used to retrieve the wanted
            %phase information.
            
            points = AOdOTF.calibrate_pupil_mask;
            [sizex,sizey] = size(AOdOTF.OTF0);
            xx = linspace(1,sizex,sizex);
            yy = linspace(1,sizey,sizey);
            [X,Y] = meshgrid(xx,yy);
            
            RA = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
            A = RA<=AOdOTF.pupil_radius;
            
            RB = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
            B = RB<=AOdOTF.pupil_radius;
            
            mask = A+B;
            mask(mask>0)=1;
            AOdOTF.plotMask = mask;
            AOdOTF.Mask = A & ~B;  
            AOdOTF.resize_mask_to_Pupil;
        end %mkMask
        
        
        function AOdOTF = centerpsfs(AOdOTF)
            if isempty(AOdOTF.psf_center)
                fprintf('Pick the center of the PSF\n');
                AOdOTF.plotPSFframe;
                p1 = pickPoint;
                AOdOTF.psf_center = p1;
            else
                p1 = AOdOTF.psf_center;
            end
            
            if isequal(size(AOdOTF.PSF0),size(AOdOTF.PSF1))
                AOdOTF.PSF0 = circshift(AOdOTF.PSF0,1-p1);
                AOdOTF.PSF1 = circshift(AOdOTF.PSF1,1-p1);
            end
        end %centerpsfs
        
        
        function AOdOTF = decenterpsfs(AOdOTF)
            p1 = AOdOTF.psf_center;
            if isequal(size(AOdOTF.PSF0),size(AOdOTF.PSF1))
                AOdOTF.PSF0 = circshift(AOdOTF.PSF0,-1+p1);
                AOdOTF.PSF1 = circshift(AOdOTF.PSF1,-1+p1);
            end
        end %decenterpsfs
        
        
        function AOdOTF = useData(AOdOTF)
            if AOdOTF.usedata == true
                AOdOTF.centerpsfs;
                AOdOTF.mkOTF;
                AOdOTF.mkdOTF;
                AOdOTF.show;
                AOdOTF.storedOTF;
            else
                warning('This command is to be used for input images of PSFs');
            end
        end %useData
        
        function AOdOTF = useData2(AOdOTF,Num_Folders,Num_files_per_folder,varargin)
            if AOdOTF.usedata == true
%                 testbedPSFs = BatchRead(Num_Folders,Num_files_per_folder, false, varargin{1,2});
%                 img_Finger = testbedPSFs{1};
%                 img_Finger = AddImages(img_Finger);
%                 
%                 imagesc(img_Finger);
%                 sqar;
%                 centerpoint = pickPoint(1);
%                 
%                 img_Finger = img_Finger(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
%                 img_Finger = img_Finger(1:end-1,1:end-1);
%                 img_No_Finger = testbedPSFs{2};
%                 img_No_Finger = AddImages(img_No_Finger);
%                 img_No_Finger = img_No_Finger(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
%                 img_No_Finger = img_No_Finger(1:end-1,1:end-1);
%                 
%                 AOdOTF.PSF0 = img_No_Finger;
%                 AOdOTF.PSF1 = img_Finger;
%                 
%                 imagesc(img_Finger);
%                 sqar;
%                 centerpoint = pickPoint(1);
%                 
%                 img_Finger = circshift(img_Finger,1-centerpoint);
%                 img_No_Finger = circshift(img_No_Finger,1-centerpoint);
                
                testbedPSFs = BatchRead(Num_Folders,Num_files_per_folder, false, varargin{1,2});
                img_Finger = testbedPSFs{1};
                img_Finger = AddImages(img_Finger);
                img_No_Finger = testbedPSFs{2};
                img_No_Finger = AddImages(img_No_Finger);
                
                
                imagesc(img_Finger);
                sqar;
                centerpoint1 = pickPoint(1);
                
                
                img_Finger = circshift(img_Finger,1-centerpoint1);
                img_Finger = fftshift(img_Finger);
                img_No_Finger = circshift(img_No_Finger,1-centerpoint1); 
                img_No_Finger = fftshift(img_No_Finger);
                
                
                
                imagesc(img_Finger);
                sqar;
                centerpoint = pickPoint(1);
                                
                img_Finger = img_Finger(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
                img_Finger = img_Finger(1:end-1,1:end-1);
                img_No_Finger = img_No_Finger(centerpoint(1) - 128:centerpoint(1) + 128,centerpoint(2) - 128:centerpoint(2) + 128);
                img_No_Finger = img_No_Finger(1:end-1,1:end-1);
                
                AOdOTF.PSF0 = img_No_Finger;
                AOdOTF.PSF1 = img_Finger;
                
                centerpoint2 = [129,129];
                img_Finger = circshift(img_Finger,1-centerpoint2);
                img_No_Finger = circshift(img_No_Finger,1-centerpoint2);
                
                
                
                
                OTF_Finger = fftshift(fft2(img_Finger));
                OTF_Finger(centerpoint2(1),centerpoint2(2)) = 0;
                OTF_No_Finger = fftshift(fft2(img_No_Finger));
                OTF_No_Finger(centerpoint2(1),centerpoint2(2)) = 0;
                
                AOdOTF.OTF0 = OTF_No_Finger;
                AOdOTF.OTF1 = OTF_Finger;
                
                AOdOTF.dOTF = AOdOTF.OTF1 - AOdOTF.OTF0;
                AOdOTF.Phase = angle(AOdOTF.dOTF);
            end
        end %useData2
        
        function scanNumericalDefocus(AOdOTF,alpha_min,alpha_max,numpoints)
            fig1 = figure(1);
            axis off;
            sqar;
            
            input('Press Enter when Ready');
            sqar;
            winsize = get(fig1,'Position');
            winsize(1:2) = [0,0];
            A = moviein(numpoints,fig1,winsize);
            set(fig1,'NextPlot','replacechildren');
            
            imagesc(AOdOTF.Phase);
%             sqar;
            axis off;
            pt = pickPoint(1);
            [X,Y] = mkImageCoords(AOdOTF.dOTF,1,pt);
            R = sqrt(X.^2 + Y.^2);
            counter = 1;
            for alpha = linspace(alpha_min,alpha_max,numpoints)
                defocus = exp(alpha*1i.*R);
                newphase = angle(defocus .* AOdOTF.dOTF);
                imagesc(newphase);
                axis off;
                title(sprintf('alpha = %0.4f ',alpha));
                drawnow;
                A(:,counter) = getframe(fig1,winsize);
                counter = counter + 1;
            end
            
            save scanDefocus.mat A winsize;
        end
        
        function scanBinning(AOdOTF,bin_min,bin_max)
            figure;
            for binning = bin_min:bin_max
                dOTF_binned = downsampleCCD(AOdOTF.dOTF,binning,binning);
                plotComplex(dOTF_binned,2);
                drawnow;
                sqar;
                title(sprintf('dOTF Signal with %d x %d binning',binning, binning));
                input('continue?');
            end
        end%scanBinning
        
        
        function [Tipdiff,Tiltdiff] = calculateTT(AOdOTF,pixelshift)
            g = AOdOTF.dOTF;
            if length(pixelshift) == 1
                pixelshift = [pixelshift,pixelshift];
            end
            tipleft = circshift(g,[0,pixelshift(2)]);
            tipright = conj(circshift(g,[0,-pixelshift(2)]));
            tiltup = circshift(g,[-pixelshift(1),0]);
            tiltdown = conj(circshift(g,[pixelshift(1),0]));
            
            Tipdiff = tipleft .* tipright;
            Tiltdiff = tiltup .* tiltdown;
        end %calculateTT
        
        
        function AOdOTF = storePSFimages(AOdOTF,image1,image2)
            AOdOTF.storedPSF0{AOdOTF.psfimage_counter} = image1;
            AOdOTF.storedPSF1{AOdOTF.psfimage_counter} = image2;
            AOdOTF.psfimage_counter = AOdOTF.psfimage_counter + 1;
        end %storePSFimages
        
        
        function AOdOTF = storedOTF(AOdOTF)
            AOdOTF.storeddOTF{AOdOTF.dOTFimage_counter} = AOdOTF.dOTF;
            AOdOTF.dOTFimage_counter = AOdOTF.dOTFimage_counter + 1;
        end %storedOTF
        
        
        function AOdOTF = clearPSFOTF(AOdOTF)
            AOdOTF.PSF0 = [];
            AOdOTF.PSF1 = [];
            AOdOTF.OTF0 = [];
            AOdOTF.OTF1 = [];
        end %clearPSFOTF
        
        
        function AOdOTF = cleardOTF(AOdOTF)
            %clears some of the properties
            AOdOTF.PSF0 = [];
            AOdOTF.PSF1 = [];
            AOdOTF.OTF0 = [];
            AOdOTF.OTF1 = [];
            AOdOTF.dOTF = [];
        end %cleardOTF
        
        
        function AOdOTF = show(AOdOTF)
            % A simple plotting tool that plots the stored PSFs, OTFs, the
            % magnitude of the dOTF, and the phase of the dOTF
            
            subplot(2,2,1)
            psf0 = AOdOTF.PSF0;
            psf1 = AOdOTF.PSF1;
            imagesc([(psf0).^0.5,(psf1).^0.5]);
            axis xy;
            axis off;
            bigtitle('PSF',12);
            daspect([1,1,1]);
            
            subplot(2,2,2)
            otf0 = AOdOTF.OTF0;
            otf1 = AOdOTF.OTF1;
            imagesc([abs((otf0).^0.5),abs((otf1).^0.5)]);
            axis xy;
            axis off;
            bigtitle('OTF',12);
            daspect([1,1,1]);
            
            subplot(2,2,3)
            imagesc(abs(AOdOTF.dOTF).^0.5);
            bigtitle('dOTF Magnitude',12);
            daspect([1 1 1]);
            axis xy;
            axis off;

            subplot(2,2,4)
            imagesc(AOdOTF.Phase);
            bigtitle('dOTF Phase',12);
            daspect([1 1 1]);
            axis xy;
            axis off;
            colormap(gray)
            drawnow;
        end %show
        
    end %methods
    
end

classdef AOScreen < AOGrid
	% AOScreen class.
	%
	% This is essentially a simplified AOPhaseScreen class, but dealing
	% only in wavefronts measured in meters.  There are methods to convert
	% the result to phase, but you have to define a reference lambda.
	%
	% Written by: Johanan L. Codona, Steward Observatory: CAAO
	% 20090417 JLCodona.  AOSim2.
	
	% Public properties
	properties(Access='public')
		altitude = 0.;	% Default is on the ground.
		lambdaRef = AOField.VBAND;
	
		mirror = false;  % This is like a height doubler.
        conjugate = false; % This changes the sign of the phase when multiplying an AOField.
		touched = true;
        radius = 1; % This is for Zernike reference.
        
        fixLF = false; % This is to patch the FFT generation using fractal scaling.
	end
	
	% Private
	properties(SetAccess='private', GetAccess='public')
		thickness;
		r0 = 0.15;
		Cn2;
        L0 = 30;
	end
	
	%% Methods
	methods
		% Constructors
		function PS = AOScreen(varargin)
			% PS = AOScreen(size_defn,[r0],[reference wavelenth])
            % This is the constructor.  It can take several arguments.
            % size_defn can be another object, the AOScreen will be built
            % to match it.
            % It can be a single number whilch will be the size of a square
            % grid with the default pixel size of 4cm.  This can all be
            % changed later and the screen rerendered.
            % The size_defn can be [n1 n2] and a non-square array will be
            % rendered.
            % r0 (optional, defaults to 15cm at VBAND) depends on a ref lambda, so you 
            % may want to define that also.
            
			% (nxy,r0,lambda)
			PS = PS@AOGrid(varargin{1});

			if(length(varargin)>1)
				PS.r0 = varargin{2};
			end
			
			if(length(varargin)>2)
				PS.lambdaRef = varargin{3};
			end

			PS.thickness = 100;
			PS.Cn2 = PS.r0^(-5/3)/0.423/(2*pi/PS.lambdaRef)^2/PS.thickness;
			
			PS.nanmap = 1; % Make clear beyond the turbulence.
		end
		
		function PS = setCn2(PS,cn2,thick)
			if(nargin<3)
				thick = 1;
			end
			
			PS.thickness = thick;
			PS.Cn2 = cn2;
			PS.r0 = (0.423*(2*pi/PS.lambdaRef)^2*cn2*PS.thickness)^(-3/5); % see Roddier
			PS.touch;
		end
		
		function PS = setR0(PS,r0,thick)
			if(nargin<3)
				thick = 1;
			end
			
			PS.thickness = thick;
			PS.r0 = r0;
			PS.Cn2 = (r0^(-5/3))/0.423/(2*pi/PS.lambdaRef)^2/PS.thickness;
			PS.touch;
        end

        function PS = setOuterScale(PS,L0)
			PS.L0 = L0;
			PS.touch;
        end        
        
		% Utilities
		
		function b = isPhase(G)
            b = true; % This is for confused programs that think this is an AOPhaseScreen.
		end
		
		function PSI = phasor(PS,lambda)
			tX(PS);
			
			if(nargin<2)
				lambda = PS.lambdaRef;
			end

			PSI = exp((2*pi*1i/lambda)*PS.grid);
		end
		
		function PHASE = phase(PS,lambda)
			tX(PS);
			
			if(nargin<2)
				lambda = PS.lambdaRef;
			end

			PHASE = (2*pi/lambda) * PS.grid;
		end
		
		function Rf = fresnel(a)
			Rf = sqrt(a.lambdaRef * a.altitude);
        end
		
        function S = touch(S)
            S.touch@AOGrid;
            S.touched = true;
        end
		
        function S = importAPP(S,FITSNAME,FRAME)
            % This is JLC specific.  Read in a PAC APP phase pattern.
            
            if(nargin<3)
                FRAME = 1;
            end
            
            HEADER = fits_read_header(FITSNAME);
            S.resize(double(HEADER.NAXIS1)); % I only support square grids for now.
            
            START = double([1,1,FRAME]);
            STOP = double([HEADER.NAXIS1,HEADER.NAXIS2,FRAME]);
            if(nargin<3)
                APP = fits_read_image(FITSNAME);
            else
                APP = fits_read_image_subset(FITSNAME,START,STOP);
            end
            % Notice the transpose.  This may cause trouble if you load
            % things from other places and don't also transpose there. I
            % will probably change this in the future.
            S.grid_ = APP'*(S.lambdaRef/2/pi); 
            %S.spacing(HEADER.XSPACING);
        end
        
        function S = addRipple(S,K,amp,phase)
        % S = addRipple(S,K,amp,phase)
            [X,Y] = S.COORDS();
            S + amp*cos(K(1)*X + K(2)*Y + phase);
            touch(S);
        end
        
        function S = addGaussian(S,CENTER,amp,width)
        % S = addGaussian(S,CENTER,amp,width)
            [X,Y] = S.COORDS();
            S + amp*exp(-((X-CENTER(1)).^2 + (Y-CENTER(2)).^2)/width^2);
            touch(S);
        end

        function S = addGDelta(S,CENTER,amp,width,xy)
            [X,Y] = S.COORDS();
            if(nargin<5)
                error('need all args specified: addGDelta(S,CENTER,amp,width,xy)');
            end

            X = X - CENTER(1);
            Y = Y - CENTER(2);
            
            if(xy==1) % x is 1.  I know.  Confusing.
                S + amp*(X.*exp(-(X.^2 + Y.^2)/width^2));
            else
                S + amp*(Y.*exp(-(X.^2 + Y.^2)/width^2));
            end
            
            touch(S);
        end

        function S = addZernike(S,n,m,amp,D,cenx,ceny)
        % AOScreen S = S.addZernike(n,m,amp,D,cenx,ceny)
            if(nargin>4)
                S.radius = D/2;
            end
            if(nargin>5)
                Xcen=cenx;    %allows segment aberrations
                Ycen=ceny;
            else	
				Xcen=0;
                Ycen=0;
			end
            
            [x,y] = S.COORDS(); % lowercase on purpose to match Zernike vars.
            %fprintf('Xcen=%d   YCen=%d \n',Xcen,Ycen);
            x = (x-Ycen)/S.radius;
            y = (y-Xcen)/S.radius;
            zern = ZernikeStringR(n,m);
            S + amp*eval(zern);
            %touch(S);
        end

		function S = addDiskHarmonic(S,n,m,amp,D,cenx,ceny)
            if(nargin>5)
                Xcen=cenx;  %allows segment aberrations
                Ycen=ceny;
            else	
				Xcen=0;
                Ycen=0;
            end
           
            if(nargin<5)
                r = S.radius;
			else
				r = D/2;
            end
            
			if(nargin<4)
                amp = S.lambdaRef/2/pi;
			end
			
            [X,Y] = S.COORDS(); 
            X = (X-Xcen)/r;
            Y = (Y-Ycen)/r;
			R = sqrt(X.^2+Y.^2);
			TH = atan2(Y,X);
            %fprintf('Xcen=%d   YCen=%d \n',Xcen,Ycen);
			
            DH = dh_dh(n,m,R,TH);
            S + amp*DH;
            touch(S);
        end
		
        function grid = LPF(S,scale)
            % grid = SCREEN.LPF(S)
            % This function returns a Gaussian smoothed version
            % of the SCREEN's displacement grid. 
            % The SCREEN itself is not altered.
            % A cheesy AO model is 
            %
            % z0 = SCREEN.make.grid;
            % z1 = SCREEN.LPF(l_actuator);
            % SCREEN.grid(1.414*(z0-z1));
            % FIELD.planewave*SCREEN*APERTURE;
            % PSF = FIELD.mkPSF(FOV,dFOV);
            % ta dah! (note that this is only a fitting error model)
            % To add lag error, shift z1 by some pixels of length
            % wind*lagtime.
        
            grid = S.grid;
            N = round(2*scale/S.dx);
            N = N + mod(N+1,2);
            filter1d = chebwin(N);
            FILTER = filter1d*filter1d';
            FILTER = FILTER / sum(FILTER(:));
            grid = conv2(grid,FILTER,'same');
            
        end
        
        function grid = LPF_(S,scale,transition)
            % grid = SCREEN.LPF(S)
            % This function returns a Gaussian smoothed version
            % of the SCREEN's displacement grid. 
            % The SCREEN itself is not altered.
            % A cheesy AO model is 
            %
            % z0 = SCREEN.make.grid;
            % z1 = SCREEN.LPF(l_actuator);
            % SCREEN.grid(1.414*(z0-z1));
            % FIELD.planewave*SCREEN*APERTURE;
            % PSF = FIELD.mkPSF(FOV,dFOV);
            % ta dah! (note that this is only a fitting error model)
        
            KMAX = 2*pi/scale;
            
            if(nargin<3)
                EXTENT = S.extent;
                %transition = 2*pi/EXTENT(1)/25
                transition = KMAX/20
            end
            
            grid = S.grid;
            [KX,KY] = S.KCOORDS;
            KR = sqrt(KX.^2+KY.^2);
            FILTER = (smoothedge(2*pi/scale-KR,transition));
            FILTER = FILTER / sum(FILTER(:));
            FILTER = fftshift(FILTER);
            
            grid = real(ifft2(fft2(S.grid).*FILTER));
            S.grid(ifft2(fft2(S.grid).*FILTER));
            %grid = conv2(grid,FILTER,'same');
            
        end
        
		function [dPhase_meanSquare,s,...
                dPhase_meanSquareSigma,...
                dPhase_] ...
                = SFestimate(PS,APER,Npoints,dspacing,lambda)
            % [dPhase_meanSquare,s,dPhase_meanSquareSigma,dPhase_] ...
            %   = SFestimate(PS,APER,Npoints,dspacing,lambda)
            % Estimate the AOScreen structure function by brute force.
            % This method Starts with Npoints in the PS grid and then
            % filters them based on being inside the pupil.
            % The number of points is reduced from the value set, and
            % pairwise ombined, so be aware of huge tasks being
            % inadvertently set.

            if(nargin<5)
                lambda = PS.lambdaRef;
            end
            
            k = 2*pi/lambda;
            
            if(nargin<4)
                ds = min(PS.spacing);
            else
                ds = dspacing;
            end
			BBox = APER.BBox;
			xymin = BBox(1,:);
			xymax = BBox(2,:);
			midpnt = mean(BBox);
            D = mean(xymax-xymin);
			
			POINTS = rand(Npoints,2);
			POINTS(:,1) = POINTS(:,1)*(xymax(1)-xymin(1))+xymin(1);
			POINTS(:,2) = POINTS(:,2)*(xymax(2)-xymin(2))+xymin(2);
		
            valid=APER.isInside(POINTS);
			
            POINTS = POINTS(valid,:);
            PHASE = k*PS.interpGrid(POINTS(:,1),POINTS(:,2));

            INDEX = 1:size(POINTS,1);

            [N,M] = meshgrid(INDEX,INDEX);
            N = N(:);
            M = M(:);
            
            dPHASE = PHASE(N) - PHASE(M);            
            dPOINTS = POINTS(N,:) - POINTS(M,:);

            S = sqrt(sum(dPOINTS.^2,2));
            sindex = int32(S/ds)+1;
    		Nbins= max(sindex);
            
            s = (double(1:Nbins)-0.5)*ds;
            
            dPhase_ = zeros(1,Nbins);
            dPhase_meanSquare = zeros(1,Nbins);
            dPhase_meanSquareSigma = zeros(1,Nbins);

            for(n=1:Nbins)
                data = dPHASE(n==sindex);
                dPhase_(n) = mean(data);
                dPhase_meanSquare(n) = mean(data.^2);
                dPhase_meanSquareSigma(n) = std(data.^2);
            end
            
            %semilogy
            loglog(s,dPhase_meanSquare,'.',...
                s,dPhase_meanSquare+dPhase_meanSquareSigma,'r--',...
                s,dPhase_meanSquare-dPhase_meanSquareSigma,'r--');
            hold on;
            plot(s([1 end]),[1 1]*6.88,'k--');
            plot([1 1]*PS.r0,[0.1 100],'k-');            
            plot(s,6.88*(s/PS.r0).^(5/3),'g--');
            hold off;
            drawnow;
        end
    
        function SF = SFtheo(PS,x)
            SF = 6.88 * (x/PS.r0).^(5/3);
        end
        
        
    end % public methods.
end

%% Make Testbed Elements
% Only run in testbed_sim, needs flags and such set there


%% ************************************************************************
%                      Construct System Elements
%**************************************************************************
if RunSIM == true
    load alright.mat;
    John = audioplayer(y,Fs);
%     play(John);
    
    %% Pupil Mask
    PUPIL_DEFN = [
        0 0 D         1 aa 0 0 0 0 0
%         0 0 secondary 0 aa/2 0 0 0 0 0
%         0 0 spider   -2 aa 4 0 D/1.9 0 0
        ];
    
    A = AOSegment;
    A.spacing(SPACING);
    A.name = 'UAWFC Pupil Mask';
    A.pupils = PUPIL_DEFN;
    A.make;
    [x,y] = A.coords;
    fprintf('\nPupil Mask Constructed\n');
    
    if system_verbose == true
        figure(1)
        A.show;
        colormap(gray);
        drawnow;
        pause(1);
    end
end

%% IrisAO DM
if RunSIM == true
    if IrisAO_on == true
        segpitch = 606e-6; %leave this alone
        magnification = 1; %leave this alone too
        FoV = FoV_withIrisAO;
        % Make the Mirror
        if Scalloped_Field == true
            [DM1,F_scal] = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
            F_scal.grid(padarray(F_scal.grid,[ceil(271/2),ceil(207/2)]));
        else
            DM1 = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
        end
        DM1.lambdaRef = lambda;
        
        % clc; %clears the command window of the text generated from adding segments to the AOAperture object
        fprintf('\nIrisAO Mirror Constructed\n');
        pupil_blocker = [1.788e-3,0.3212e-3,0.1e-3];
    else
        DM1 = 1;
        FoV = FoV_withoutIrisAO;
        Scalloped_Field = false;
        fprintf('\nUsing a Flat instead of the IrisAO\n');
        pupil_blocker = [3.28e-3,0.3212e-3,0.1e-3];
    end
end

%% Set the Initial Piston, Tip, Tilt of IrisAO Mirror
% Flatten the IrisAO
PTTpos = zeros(37,3);

% Set a Random Mirror Shape
% PTTpos = horzcat(randn(37,1)*10^-6,randn(37,1)*10^-3,randn(37,1)*10^-3);

% Use a Zernike on the mirror....magnification must be equal to 1 for this!
% Zernike_Number = [2,3];
% Zernike_Coefficient_waves = randn(2,1);
% PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);

if RunTESTBED == true
    fprintf('*******************************\nSending PTTpos to IrisAO Mirror\n*******************************\n\n');
    %     tempdir = pwd;
    %     cd /home/lab/Desktop/Shared_Stuff
    %     save('PTTpos','PTTpos');
    %     checkmirrorupdated = false;
    %     while(checkmirrorupdated == false)
    %         CMD_FILES = dir('PTTpos.mat');
    %         if(~isempty(CMD_FILES))
    %             pause(0.1);
    %         else
    %             checkmirrorupdated = true;
    %         end
    %     end
    %     cd(tempdir);
    %     clear tempdir;
end

if RunSIM == true
    if IrisAO_on == true
        % Load in Mapping Data
        load('IrisAO_SegMap.mat');
        PTT = zeros(37,3);
        
        % Map the PTT matrix from hardware to software order
        for ii = 1:37
            mapped_segment = IrisAO_SegMap(ii);
            PTT(ii,1:3) = PTTpos(mapped_segment,:);
        end
        
        % Send to DM Model
        DM1.PTT(PTT);
        DM1.touch;
        DM1.render;
        
        if system_verbose == true
            figure(1);
            DM1.show;
            colormap(gray);
            drawnow;
            pause(1);
        end
    end
end

%% Boston MircroMachines DM
%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
STROKE = 1.5e-6;
Pitch = 300e-6;
radius_BMC = (4.65)*10^-3;

if RunSIM == true
    if BMC_on == true
        %Coordinate System
        xmin = -radius_BMC;
        xmax = radius_BMC;
        BMC_x = (xmin:SPACING:xmax);
        ymin = -radius_BMC;
        ymax = radius_BMC;
        BMC_y = (ymin:SPACING:ymax);
        [BMC_X,BMC_Y] = meshgrid(BMC_x,BMC_y);
        
        BMC_pupil = ones(size(BMC_X));
        % BMC_pupil = padarray(BMC_pupil,[250,250],1,'both'); %put 250 pixels on both sides outside of active area
        
        %Construct the Pupil
        Seg = AOSegment(length(BMC_pupil));
        Seg.spacing(SPACING);
        Seg.name = 'BMC Pupil';
        Seg.grid(BMC_pupil);
        
        %Make it an AOAperture Class
        A_BMC = AOAperture;
        A_BMC.spacing(SPACING);
        A_BMC.name = 'BMC Aperture';
        A_BMC.addSegment(Seg);
        
        %Make an AODM
        DM2 = AODM(A_BMC);
        [X,Y] = DM2.COORDS;
        
        % Create Actuator Locations
        actuator_x = xmin:Pitch:xmax;
        actuator_x = actuator_x - mean(actuator_x);
        actuator_y = ymin:Pitch:ymax;
        actuator_y = actuator_y - mean(actuator_y);
        [ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
        %     ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
        %     ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
        %     ACTUATOR_X(ACTUATOR_X==0) = [];
        %     ACTUATOR_Y(ACTUATOR_Y==0) = [];
        BMC_ACTS = [ACTUATOR_X(:),ACTUATOR_Y(:)];
        nacts = length(BMC_ACTS(:,1));
        
        % Add Actuators
        DM2.addActs(BMC_ACTS,1,A_BMC);
        
        % Turn Off Actuators that Aren't Illuminated
        RHO = zeros(nacts,1);
        for ii = 1:nacts
            RHO(ii) = sqrt(BMC_ACTS(ii,1)^2 + BMC_ACTS(ii,2)^2);
            if RHO(ii) > D/2 + D/18
                if RHO(ii) < ((D/2) + (D/18))
                    DM2.actuators(ii,5) = 2;
                else
                    DM2.actuators(ii,5) = 0;
                end
            elseif RHO(ii) < secondary/2.1
                %             DM2.actuators(ii,5) = 0;
            end
        end
        
        % Get List of Which Actuators are Being Used
        DM2.setOnActs;
        
        %Turn Off Actuators so the Program Knows they are off
        DM2.disableActuators(DM2.OffActs);
        
        % Set the Convex Hull Boundary Conditions
        %     DM2.defineBC(D/2 + D/18,4,'circle');
        DM2.defineBC(D/2,4,'circle');
        
        [x_DM,y_DM] = DM2.coords;
        %% Set the Initial Piston Values of the BMC DM
        DM2.flatten;
        %     DM2.actuators(343,3) = 1e-6;
        %     DM2.actuators(OnActs,3) = ((randn(length(OnActs),1).*10^-6));
        DM2.removeMean;
        
        if RunTESTBED == true
            DM_Pistons = reshape(DM2.actuators(:,3),[32,32]);
            
            tempdir = pwd;
            cd C:\Users\atrod_000\Documents\GitHub\UAWFC\fits_files;
            fitswrite(DM_Pistons,'DM_Pistons.fits');
            img = fitsread('DM_Pistons.fits');
            figure(5);
            imagesc(img);
            input('Press Enter');
            close
            cd(tempdir);
            
            % DO STUFF TO SEND TO MIRROR
            
            
        end
        
        if system_verbose == true
            DM2.show;
            colormap(gray);
            title('BMC Pupil');
            DM2.plotActuators;
            pause(1);
            DM2.plotRegions;
            pause(1);
        end
        fprintf('\nBMC DM Constructed\n');
    else
        DM2 = 1;
        [X,Y] = A.COORDS;
        fprintf('\nUsing a Flat instead of the BMC\n');
    end
    
end
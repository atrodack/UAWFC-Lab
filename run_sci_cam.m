%run_cam.m
%________________________________________________________________________
%Contact: K.Miller [millerk2@email.arizona.edu]               August 2014
%Contact: A.Rodack [atrodack@email.arizona.edu]
%Contact: J.Knight [jknight@optics.arizona.edu]
%
%Opens imaqtool, identifies all cameras and sets camera settings. Can run
%a live video preview or take data.  For the latter, must be in
%/home/genius2/Desktop directory.  Allows for taking dark and background
%images for improved final image quality. To use previous dark and
%background images, they must be stored in the Data folder, NOT in
%an individual batch folder.
%_________________________________________________________________________
clear all;
clc;
close all;

%% Create dOTF Object
cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
directory_dOTFclass = pwd;
nxy = 1;  %This shouldn't matter, but we will find out when this code is tested
dOTF = AOdOTF(nxy);
counterloop = 1;
%% Open imaqtool
imaqtool

%% Camera Acquistion and Settings
%Recognizes and labels FPM and Science Image and sources
scim = videoinput('dcam',2,'Y16_640x480'); %updated adaptorname to match matlab adaptorname from imaqhwinfo on 2/4/2015 by Justin
% srcscim=getselectedsource(scim);
% 
% %Sets camera settings
% %src.GammaMode = 'manual';
% %src.Gamma = 0;
% scim.FramesPerTrigger=1;
% srcscim.BusSpeed='S400';
% srcscim.FrameRate='30';

disp('Run preview   [1]');
disp('Take images   [2]');
select=input(':');
restartcam = 1;
p2 = [];
%% loop the program
while(restartcam == 1)
    cd /home/lab/Desktop
    %% Clear All Existing variables, previews and figures
    clearvars -except dOTF counter_ counterloop restartcam select scim p2 p1 p3 p7 p8 p9 points
    close all
    closepreview
    clc
srcscim=getselectedsource(scim);

%Sets camera settings
%src.GammaMode = 'manual';
%src.Gamma = 0;
scim.FramesPerTrigger=1;
srcscim.BusSpeed='S400';
srcscim.FrameRate='30'; 

    
    %% Camera Preview and Data Acquisition
    switch(select)
        
        case 1
            %% Preview Case
            preview(scim);
            %         preview(fpm);
            %Choose shutter speed
            disp(' ')
            disp('Science image shutter speed (1-531)')
            disp('(Suggested value for viewing: 5)')
            disp('(Suggested value: 10 - 15 without ND in, 70 with)');
            sSpeed=input(': ');
            srcscim.Shutter=sSpeed;
            restartcam = 0;
        case 2
            %% Data Acquistion Case
            preview(scim);
            %         preview(fpm);
%             if counterloop == 1
                %Choose shutter speed
                disp(' ')
                disp('Science image shutter speed (1-531)')
                disp('(Suggested value: 10 - 15)');
                sSpeed=input(': ');
                srcscim.Shutter=sSpeed;
%             end
            
            %Sets number of frames to be taken
            disp(' ')
            nframes=input('Number of frames: ');
            
            %Changes into Data folder
            cd Data
            
            %Create variable of all files/folders in Data Directory
            listing=ls('/home/lab/Desktop/Data');
            %listing=ls('/home/genius2/Desktop/Data');
            
            %% Ask if User wants to use previously taken darks and backgrounds
            fprintf('Do you want to use Previously Created Dark and Background Images?\n');
            bkgrnddarkload = input('[1] for Yes, [0] for No: ');
            
            %Take in folder that has saved images, create string of its
            %location, add it to path, check for image existance, remove from
            %path
            if bkgrnddarkload == 1
                disp('--------------------------------------------------------------')
                disp(listing)
                disp('--------------------------------------------------------------')
                disp(' ')
                dn=input('Copy Folder Name Here that has Dark and Background images: ','s');
                dirnamebkgrnddark=(['/home/lab/Desktop/Data',dn]);
                %dirnamebkgrnddark=(['/home/genius2/Desktop/Data/',dn]);
                addpath(dirnamebkgrnddark);
                %Checks for previously loaded dark and background images
                cSb=exist('SCIbkgd.fits','file');
                rmpath(dirnamebkgrnddark);
                %Set directory name to the Data folder.  The exist check will
                %return 0. This ensures the code works correctly further down
            else
                dirnamebkgrnddark = '/home/lab/Desktop/Data';
                cSb=exist('SCIbkgd.fits','file');
            end
            
            %% Create folder for each run
            disp(' ')
            %Ask user for data run number
            B=input('Batch Number: ');
            %Find date from computer clock
            t=datevec(datestr(clock,0));
            
            disp(' ')
            %Ask user if a filter is in place
            disp('Imaging with filter?: ')
            disp('BP 40nm       [1]')
            disp('BP 10nm       [2]')
            disp('No filter     [3]')
            filter=input(': ');
            if filter==1
                filterselection=sprintf('bp40');
            end
            if filter==2
                filterselection=sprintf('bp10');
            end
            if filter==3
                filterselection=sprintf('nofilter');
            end
            
            disp(' ')
            %Ask user to input additional title identifier
            disp('Additional title identifier?: ')
            disp('Yes           [1]')
            disp('No            [0]')
            extratitle=input(':');
            
            %Creat folder with above user inputs
            switch(extratitle)
                case 1
                    %% Additional Title Identification
                    %get identifier from user
                    et=input('Title identifier: ','s');
                    %make the folder
                    mkdir([num2str(t(1)),num2str(t(2)),num2str(t(3)),'_Batch',num2str(B),'_',filterselection,'_',et]);
                    %store folder name to variable
                    dirnamebatchfolder = [num2str(t(1)),num2str(t(2)),num2str(t(3)),'_Batch',num2str(B),'_',filterselection,'_',et];
                    %Change into folder that has Previously Created Dark and Background IMages
                    cd(dirnamebkgrnddark)
                    if cSb ==2
                        copyfile('SCIbkgd.fits',['/home/lab/Desktop/Data/',dirnamebatchfolder])
                    end
                    %Change back into created folder for current run
                    cd (['/home/lab/Desktop/Data/',dirnamebatchfolder])
                    
                case 0
                    %% No Additional Title Identification
                    %make the folder
                    mkdir([num2str(t(1)),num2str(t(2)),num2str(t(3)),'_Batch',num2str(B),'_',filterselection]);
                    %store folder name to variable
                    dirnamebatchfolder = [num2str(t(1)),num2str(t(2)),num2str(t(3)),'_Batch',num2str(B),'_',filterselection];
                    %Change into folder that has Previously Created Dark and Background IMages
                    cd(dirnamebkgrnddark)
                    %Copy them into the created folder if the above code says they exist
                    if cSb ==2
                        copyfile('SCIbkgd.fits',['/home/lab/Desktop/Data/',dirnamebatchfolder])
                    end
                    %Change back into created folder for current run
                    cd (['/home/lab/Desktop/Data/',dirnamebatchfolder])
                    
            end
            %% Take or Load in Dark and Background Images
            %If no was selected earlier
            if bkgrnddarkload ~= 1
                disp(' ')
                disp('Take background frame?')
                disp('Yes           [1]')
                disp('No            [0]')
                bg=input(':');
                if bg==1
                    disp('Press ENTER when ready: ')
                    pause
                    start(scim);
                    SCIbkgd=getdata(scim,1);
                    SCIbkgd=int16(SCIbkgd/1024);
                    fitswrite(SCIbkgd,'SCIbkgd.fits');
                else
                    FPMbkgd=0;
                    SCIbkgd=0;
                end
            else
                disp(' ')
                disp('Use preloaded background?')
                disp('Yes           [1]')
                disp('No            [0]')
                new=input(':');
                if new==1
                    bg = 1;
                    if cFb == 0
                        disp('WARNING: No FPM background image saved')
                        FPMbkgd=0;
                    else
                        
                    end
                    if cSb == 0
                        disp('WARNING: No science background image saved')
                        SCIbkgd=0;
                    else
                        %                     SCIbkgd=uint8(fitsread('SCIbkgd.fits'));
                        SCIbkgd=fitsread('SCIbkgd.fits');
                    end
                    
                else
                    FPMbkgd=0;
                    SCIbkgd=0;
                    bg = 0;
                end
            end
            
            disp(' ')
            %% Imaging
            disp('Press ENTER to begin imaging: ')
            pause
            counter=0;
            for i=1:nframes
                counter=counter+1;
                disp(counter)
                %Opens video display and begins data acquisition
                start(scim);
                
                %Stores image, timestamp, and camera settings
                [scimage, scitime, scimetadata]=getdata(scim,1);
                scimage=int16(scimage/1024);
                
                SCIcorrected=scimage-SCIbkgd;
                
                %Opens metadata files
                scimdata=struct2cell(scimetadata);
                ts=scimdata{1};
                framenum=scimdata{2};
                relativeframe=scimdata{3};
                triggerindex=scimdata{4};
                
                %Writes image file and names with date and timestamp
                %Order for image naming: YEAR MONTH DAY HOUR MINUTE
                fitswrite(scimage,['RAW_scienceIM_frame_',num2str(i),'.fits'])
                psfimage = fitsread(['RAW_scienceIM_frame_',num2str(i),'.fits']);
                
                cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
                if counterloop == 1
                    dOTF.PSF0 = psfimage;
                else
                    if mod(counter_,2) == 1
                        dOTF.PSF0 = psfimage;
                    elseif mod(counter_,2) == 0
                        dOTF.PSF1 = psfimage;
                    end
                end
                cd (['/home/lab/Desktop/Data/',dirnamebatchfolder])
                
                
                if bg==1
                    fitswrite(SCIcorrected,['scienceIM_frame_',num2str(i),'.fits'])
                    psfimagecorrected = fitsread(['scienceIM_frame_',num2str(i),'.fits']);
                    
                    cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
                    if counterloop == 1
                        dOTF.PSF0 = psfimagecorrected;
                    else
                        if mod(counter_,2) == 1
                            dOTF.PSF0 = psfimagecorrected;
                        elseif mod(counter_,2) == 0
                            dOTF.PSF1 = psfimagecorrected;
                        end
                        cd (['/home/lab/Desktop/Data/',dirnamebatchfolder])
                    end
                end
                
                cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
                if ~isempty(dOTF.PSF1)
                    dOTF.storePSFimages(dOTF.PSF0,dOTF.PSF1);
                end
                
            end
            
            
            
            %% dOTF Calculations
            
            if counterloop ~= 1
                disp(' ')
                disp('Perform dOTF?')
                disp('Yes           [1]')
                disp('No            [0]')
                decision = input(':');
                
                if decision == 1
                    disp(' ')
                    disp('Use Picture 1 as non-Modified Pupil?)')
                    disp('Yes           [1]')
                    disp('No            [0]')
                    decision2 = input(':');
                    
                    if decision2 == 1
                        PSF = dOTF.storedPSF0{1};
                        if mod(counter_,2) == 1
                            dOTF.PSF1 = PSF;
                        elseif mod(counter_,2) == 0
                            dOTF.PSF0 = PSF;
                        end
                        
                        dOTF.useData;
                        dOTF.cleardOTF;
                        cd /home/lab/Desktop
                    else
                        dOTF.useData;
                        dOTF.cleardOTF;
                        cd /home/lab/Desktop
                    end
                else
                    cd /home/lab/Desktop
                end
            else
                cd /home/lab/Desktop
                decision = 0;
            end
            
            %% Send to DM
            if decision == 1
                disp(' ')
                disp(' ')
                disp('Send to DM?')
                disp('Yes           [1]')
                disp('No            [0]')
                sendtoDM = input(': ');
                if sendtoDM == 1
                    mask_radius = 5; %pixels
                    if isempty(p2)
                        clf;
                        dOTF.unwrapphase('gold');
                        imagesc(dOTF.Phase);
                        colormap(jet);
                        daspect([1,1,1]);
                        p2 = pickPoint;
                        p1 = pickPoint;
                        p3 = pickPoint;
                        p7 = pickPoint;
                        p8 = pickPoint;
                        p9 = pickPoint;
                        
                        points = cell(1,3);
                        points{1} = p2;
                        points{2} = p1;
                        points{3} = p3;
                        points{4} = p7;
                        points{5} = p8;
                        points{6} = p9;
                    end
                    
                    phase = dOTF.Phase;
                    phase = phase ./ ((2*pi)/0.6);
                    [sizey,sizex] = size(phase);
                    xx = linspace(1,sizex,sizex);
                    yy = linspace(1,sizey,sizey);
                    [X,Y] = meshgrid(xx,yy);
                    
                    R2 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
                    M2 = R2<=mask_radius;
                    seg2 = M2.*phase;
                    seg2piston = mean(seg2(:));
                    
                    R1 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
                    M1 = R1<=mask_radius;
                    seg1 = M1.*phase;
                    seg1piston = mean(seg1(:));
                    
                    R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
                    M3 = R3<=mask_radius;
                    seg3 = M3.*phase;
                    seg3piston = mean(seg3(:));
                    
                    R7 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
                    M7 = R7<=mask_radius;
                    seg7 = M7.*phase;
                    seg7piston = mean(seg7(:));
                    
                    R8 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
                    M8 = R8<=mask_radius;
                    seg8 = M8.*phase;
                    seg8piston = mean(seg8(:));
                    
                    R9 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
                    M9 = R9<=mask_radius;
                    seg9 = M9.*phase;
                    seg9piston = mean(seg9(:));
                    
                    clf;
                    imagesc(seg2+seg1+seg3+seg7+seg8+seg9);
                    axis xy;
                    daspect([1,1,1]);
                    
                    PTTpos = zeros(37,3);
                    PTTpos(2,1) = seg2piston;
                    PTTpos(1,1) = seg1piston;
                    PTTpos(3,1) = seg3piston;
                    PTTpos(7,1) = seg7piston;
                    PTTpos(8,1) = seg8piston;
                    PTTpos(9,1) = seg9piston;
                    %Convert to base unit of micron
%                     PTTpos = PTTpos*1e-6; ***Not needed because I am
%                     dividing by 0.6 not 0.6e-6

                    cd /home/lab/Desktop/Shared_Stuff
                    save('PTTpos','PTTpos');
                else
                    fprintf('Not Sending Data to DM\n');
                end
            end
                
            
            %% Restart Program
            disp(' ')
            disp(' ')
            disp('Take another image?')
            disp('Yes           [1]')
            disp('No            [0]')
            restartcam=input(': ');
            if restartcam == 1
                if counterloop == 1
                    counter_ = counterloop+1;
                    counterloop = 0;
                else
                    counter_ = counter_ + 1;
                    
                end
            else
                fprintf('Program saving dOTF_structure (if any) and closing\n');
            end
            
            
            
            
    end
    
    
    
end



if select == 2
    cd('/home/lab/Desktop/dOTF_saved_structures')
    textstring = [num2str(t(1)),num2str(t(2)),num2str(t(3)),'dOTFstructure_data_run_',num2str(B)];
    save(textstring,'dOTF')
    cd /home/lab/Desktop
end
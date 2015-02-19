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
DEBUG_FLAG = true;

%% Create dOTF Object
cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
directory_dOTFclass = pwd;
nxy = 1;  %This shouldn't matter, but we will find out when this code is tested
dOTF = AOdOTF(nxy);
counterloop = 1;
counter_ = 1;
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
total_loops_max = 4;
p2 = [];
%% loop the program
while(counter_ <= total_loops_max)
    cd /home/lab/Desktop/
    %% Clear All Existing variables, previews and figures
    clearvars -except dOTF counter_ counterloop total_loops_max restartcam select scim bg cSb SCIbkgd Datasubfolderbydate nframes dirnamebkgrnddark extratitle filterselection bkgrnddarkload B p2 p1 p3 p7 p8 p9 points DEBUG_FLAG
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
            counter_ = 1e19; %kill while loop
        case 2
            %% Data Acquistion Case
            preview(scim);
            %         preview(fpm);
                %Choose shutter speed
                disp(' ')
                disp('Science image shutter speed = 15')

                sSpeed=15;
                srcscim.Shutter=sSpeed;
                
                if counterloop == 1
                    %Sets number of frames to be taken
                    disp(' ')
                    nframes=input('Number of frames: ');
                    
                    %Changes into Data folder
                    cd /home/lab/Desktop/Data
                    
                    %Create variable of all files/folders in Data Directory
                    listing=ls('/home/lab/Desktop/Data');
                end
            
            %% Ask if User wants to use previously taken darks and backgrounds
            if counterloop == 1
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
            end
            %% Create folder for each run
            if counterloop == 1
                disp(' ')
                %Ask user for data run number
                B=input('Batch Number: ');
            end
            
            %Find date from computer clock
            t=datevec(datestr(clock,0));
            
            disp(' ')
            
            if counterloop == 1
               %make a directory with the current date/time under Data that will have
               %all of the current runs folders stored within it
               mkdir([num2str(t(1)),'_',num2str(t(2)),'_',num2str(t(3)),'_',num2str(t(4)),'_',num2str(t(5))]);
               Datasubfolderbydate = [num2str(t(1)),'_',num2str(t(2)),'_',num2str(t(3)),'_',num2str(t(4)),'_',num2str(t(5))];
            end
                cd(['/home/lab/Desktop/Data/',Datasubfolderbydate]);
                

            
            %Ask user if a filter is in place
            if counterloop == 1
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
            end
            
            disp(' ')
            extratitle=1;

            %Creat folder with above user inputs
                    %% Additional Title Identification
                    %get identifier from user
                    if mod(counter_,2) == 1
                        et = ['nofinger',num2str(counter_)];
                    else
                        et = ['test',num2str(counter_)];
                    end
                    
                    %make the folder
                    mkdir(['Batch',num2str(B),'_',filterselection,'_',et]);
                    %store folder name to variable
                    dirnamebatchfolder = ['Batch',num2str(B),'_',filterselection,'_',et];
                    %Change into folder that has Previously Created Dark and Background IMages
                    cd(dirnamebkgrnddark)
                    if cSb ==2
                        copyfile('SCIbkgd.fits',['/home/lab/Desktop/Data/',Datasubfolderbydate,'/',dirnamebatchfolder])
                    end
                    %Change back into created folder for current run
                    cd (['/home/lab/Desktop/Data/',Datasubfolderbydate,'/',dirnamebatchfolder])
            %% Take or Load in Dark and Background Images
            %If no was selected earlier
            if counterloop == 1
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
                        SCIbkgd=int16(SCIbkgd/2^4);
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
            end
            
            
            
            %% Use Segment Piston as Finger
            
            if mod(counter_,2) == 0
                segfinger = 1;
                fprintf('Using Segment 19 as Finger\n');
                lambda = 0.6; %microns
                pokeval = lambda / 4;
                PTTpos = zeros(37,3);
                PTTpos(19,1) = pokeval;
                tempdir = pwd;
                cd /home/lab/Desktop/Shared_Stuff
                save('PTTpos','PTTpos');
                cd(tempdir);
                clear tempdir;
            else
                segfinger = 0;
                fprintf('No finger case\n');
                PTTpos = zeros(37,3);
                tempdir = pwd;
                cd /home/lab/Desktop/Shared_Stuff
                save('PTTpos','PTTpos');
                cd(tempdir);
                clear tempdir;
            end
              
            %% Check if Mirror is updated
            if DEBUG_FLAG == false
                if mod(counter_,2) == 0
                    tempdir = pwd;
                    cd /home/lab/Desktop/Shared_Stuff
                    checkmirrorupdated = false;
                    while(checkmirrorupdated == false)
                        CMD_FILES = dir('PTTpos.mat');
                        if(~isempty(CMD_FILES))
                            pause(0.1);
                        else
                            checkmirrorupdated = true;
                        end
                    end
                    cd(tempdir);
                    clear tempdir;
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
                whos scimage
                scimage=int16(scimage/(2^4));
                
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
                cd (['/home/lab/Desktop/Data/',Datasubfolderbydate,'/',dirnamebatchfolder])
                
                
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
                        cd (['/home/lab/Desktop/Data/',Datasubfolderbydate,'/',dirnamebatchfolder])
                    end
                end
                
                cd('/home/lab/Desktop/Alex_AOSim2/AOSim2-rodack/AOSim2')
                if ~isempty(dOTF.PSF1)
                    dOTF.storePSFimages(dOTF.PSF0,dOTF.PSF1);
                end
                
            end
            
            
            
            %% dOTF Calculations
            
            if mod(counter_,2) == 0
                dOTF.useData;
                pause(5);
                dOTF.cleardOTF;
                cd /home/lab/Desktop/
                decision = 1;
            else
                cd /home/lab/Desktop/
                decision = 0;
            end
            
            %% Send to DM
            if decision == 1
                sendtoDM = 1;
                if sendtoDM == 1
                    mask_radius = 5; %pixels
                    if isempty(p2)
                        fprintf('Pick the centers of the useful segments\n');
                        clf;
%                         dOTF.unwrapphase('gold');
                        imagesc(abs(dOTF.storeddOTF{1}).^0.5);
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
                    dOTF.unwrapphase('gold');
                    phase = dOTF.Phase;
                    
                    %convert to OPL?
%                     phase = phase ./ ((2*pi)/0.6);

                    [sizey,sizex] = size(phase);
                    xx = linspace(1,sizex,sizex);
                    yy = linspace(1,sizey,sizey);
                    [X,Y] = meshgrid(xx,yy);
                    
                    R2 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
                    M2 = R2<=mask_radius;
                    seg2 = M2.*phase;
                    seg2Prop = seg2(seg2>0);
                    seg2piston = mean(seg2Prop);
%                     seg2piston = mean(seg2(:));
                    
                    R1 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
                    M1 = R1<=mask_radius;
                    seg1 = M1.*phase;
                    seg1Prop = seg1(seg1>0);
                    seg1piston = mean(seg1Prop);
%                     seg1piston = mean(seg1(:));
                    
                    R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
                    M3 = R3<=mask_radius;
                    seg3 = M3.*phase;
                    seg3Prop = seg3(seg3>0);
                    seg3piston = mean(seg3Prop);
%                     seg3piston = mean(seg3(:));
                    
                    R7 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
                    M7 = R7<=mask_radius;
                    seg7 = M7.*phase;
                    seg7Prop = seg7(seg7>0);
                    seg7piston = mean(seg7Prop);
%                     seg7piston = mean(seg7(:));
                    
                    R8 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
                    M8 = R8<=mask_radius;
                    seg8 = M8.*phase;
                    seg8Prop = seg8(seg8>0);
                    seg8piston = mean(seg8Prop);
%                     seg8piston = mean(seg8(:));
                    
                    R9 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
                    M9 = R9<=mask_radius;
                    seg9 = M9.*phase;
                    seg9Prop = seg9(seg9>0);
                    seg9piston = mean(seg9Prop);
%                     seg9piston = mean(seg9(:));
                    
                    
%                     clf;
%                     imagesc(seg2+seg1+seg3+seg7+seg8+seg9);
%                     daspect([1,1,1]);
                    
                    PTTpos = zeros(37,3);
                    PTTpos(2,1) = seg2piston;
                    PTTpos(1,1) = seg1piston;
                    PTTpos(3,1) = seg3piston;
                    PTTpos(7,1) = seg7piston;
                    PTTpos(8,1) = seg8piston;
                    PTTpos(9,1) = seg9piston;
                    
                    if segfinger == 1
                        PTTpos(19,1) = -pokeval;
                    end
                    PTTpos
                    input('Press Enter to Send to DM');

                    cd /home/lab/Desktop/Shared_Stuff
                    save('PTTpos','PTTpos');
                else
                    fprintf('Not Sending Data to DM\n');
                end
            end
                
                input('Press Enter to Continue');
            %% Restart Program
            
                if counterloop == 1
                    counter_ = counterloop+1;
                    counterloop = 0;
                else
                    counter_ = counter_ + 1; 
                end

            
            
            
            
    end
    
    
    
end



if select == 2
    cd('/home/lab/Desktop/dOTF_saved_structures')
    textstring = [num2str(t(1)),'_',num2str(t(2)),'_',num2str(t(3)),'_',num2str(t(4)),'_',num2str(t(5)),'_','dOTFstructure_data_run_',num2str(B)];
    save(textstring,'dOTF')
    cd /home/lab/Desktop/
end
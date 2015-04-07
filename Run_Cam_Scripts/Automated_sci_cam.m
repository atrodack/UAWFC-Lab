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
Send2DM = false;
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
total_loops_max = 1;
p1 = [];
%% loop the program
while(counter_ <= total_loops_max)
    cd /home/lab/Desktop/
    %% Clear All Existing variables, previews and figures
    clearvars -except dOTF counter_ counterloop total_loops_max restartcam Send2DM select scim bg cSb SCIbkgd Datasubfolderbydate nframes dirnamebkgrnddark extratitle filterselection bkgrnddarkload B p1 points DEBUG_FLAG
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
            counter_ = 1e19; %kill truewhile loop
        case 2
            %% Data Acquistion Case
            preview(scim);
            %         preview(fpm);
            %Choose shutter speed
            disp(' ')
            
            
            if counterloop == 1
                sSpeed= 10;
            else
                sSpeed = 10;
            end
            
            fprintf('Science image shutter speed = %g',sSpeed)
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
                fprintf('Using Segment 20 as Finger\n');
                lambda = 0.6; %microns
                pokeval = lambda / 4;
                PTTpos = zeros(37,3);
                PTTpos(20,1) = pokeval;
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
            
            %% Unpoke the pokeseg
            if segfinger == 1
                PTTpos(20,1) = -pokeval;
            end
                    
            %% Send to DM
            if decision == 1
                if Send2DM == true
                    mask_radius = 5; %pixels
                    if isempty(p1)
                        fprintf('Pick the centers of the useful segments\n');
                        clf;
                        imagesc(abs(dOTF.storeddOTF{1}).^0.5);
                        colormap(jet);
                        daspect([1,1,1]);
                        %Pick the points
                        p1 = pickPoint;
                        p2 = pickPoint;
                        p3 = pickPoint;
                        p4 = pickPoint;
                        p5 = pickPoint;
                        %p6 = pickPoint; %Segment is Dead
                        p7 = pickPoint;
                        p8 = pickPoint;
                        p9 = pickPoint;
                        p10 = pickPoint;
                        p11 = pickPoint;
                        p12 = pickPoint;
                        p13 = pickPoint;
                        p14 = pickPoint;
                        p15 = pickPoint;
                        %p16 = pickPoint; %Segment is Dead
                        p17 = pickPoint;
                        p18 = pickPoint;
                        p19 = pickPoint;
%                         p20 = pickPoint; %pokeseg
                        p21 = pickPoint;
                        p22 = pickPoint;
                        p23 = pickPoint;
                        p24 = pickPoint;
                        p25 = pickPoint;
                        p26 = pickPoint;
                        p27 = pickPoint;
                        p28 = pickPoint;
                        %p29 = pickPoint; %Segment is Dead
                        %p30 = pickPoint; %Segment is Dead
                        %p31 = pickPoint; %Segment is Dead
                        %p32 = pickPoint; %Segment is Dead
                        %p33 = pickPoint; %Segment is Dead
                        %p34 = pickPoint; %Segment is Dead
                        %p35 = pickPoint; %Segment is Dead
                        p36 = pickPoint;
                        p37 = pickPoint;
                        
                        %store the picked points
                        points = cell(1,37);
                        points{1} = p1;
                        points{2} = p2;
                        points{3} = p3;
                        points{4} = p4;
                        points{5} = p5;
                        %points{6} = p6; %Segment is Dead
                        points{7} = p7;
                        points{8} = p8;
                        points{9} = p9;
                        points{10} = p10;
                        points{11} = p11;
                        points{12} = p12;
                        points{13} = p13;
                        points{14} = p14;
                        points{15} = p15;
                        %points{16} = p16; %Segment is Dead
                        points{17} = p17;
                        points{18} = p18;
                        points{19} = p19;
%                         points{20} = p20; %pokeseg
                        points{21} = p21;
                        points{22} = p22;
                        points{23} = p23;
                        points{24} = p24;
                        points{25} = p25;
                        points{26} = p26;
                        points{27} = p27;
                        points{28} = p28;
                        %points{29} = p29; %Segment is Dead
                        %points{30} = p30; %Segment is Dead
                        %points{31} = p31; %Segment is Dead
                        %points{32} = p32; %Segment is Dead
                        %points{33} = p33; %Segment is Dead
                        %points{34} = p34; %Segment is Dead
                        %points{35} = p35; %Segment is Dead
                        points{36} = p36;
                        points{37} = p37;
                        
                    end
                    shift_point = dOTF.pupil_center;
                    [sizex,sizey] = size(dOTF.Phase);
                    center = [round(sizex/2),round(sizey/2)];
                    
                    phase_ref = dOTF.Phase(center(1),center(2));
                    
                    dOTF.unwrapphase('gold');
                    phase = dOTF.Phase - phase_ref;
                    
                    %Tip / Tilt calculation
                    
                    % convert to OPL
                    OPL = phase ./ ((2*pi)/0.55);
                    
                    [sizey,sizex] = size(OPL);
                    xx = linspace(1,sizex,sizex);
                    yy = linspace(1,sizey,sizey);
                    [X,Y] = meshgrid(xx,yy);
                    
                    R2 = sqrt((X-points{2}(2)).^2 + (Y-points{2}(1)).^2);
                    M2 = R2<=mask_radius;
                    seg2 = M2.*OPL;
                    [tip2,tilt2] = calctiptiltdOTF(M2.*phase,p2,[3,3]);
                    seg2Prop = seg2(abs(seg2)>0);
                    seg2piston = mean(seg2Prop);
                    
                    
                    R1 = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
                    M1 = R1<=mask_radius;
                    seg1 = M1.*OPL;
                    [tip1,tilt1] = calctiptiltdOTF(M1.*phase,p1,[3,3]);
                    seg1Prop = seg1(abs(seg1)>0);
                    seg1piston = mean(seg1Prop);
                    
                    
                    R3 = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
                    M3 = R3<=mask_radius;
                    seg3 = M3.*OPL;
                    [tip3,tilt3] = calctiptiltdOTF(M3.*phase,p3,[3,3]);
                    seg3Prop = seg3(abs(seg3)>0);
                    seg3piston = mean(seg3Prop);
                    
                    
                    R4 = sqrt((X-points{4}(2)).^2 + (Y-points{4}(1)).^2);
                    M4 = R4<=mask_radius;
                    seg4 = M4.*OPL;
                    [tip4,tilt4] = calctiptiltdOTF(M4.*phase,p4,[3,3]);
                    seg4Prop = seg4(abs(seg4)>0);
                    seg4piston = mean(seg4Prop);
                    
                    
                    R5 = sqrt((X-points{5}(2)).^2 + (Y-points{5}(1)).^2);
                    M5 = R5<=mask_radius;
                    seg5 = M5.*OPL;
                    [tip5,tilt5] = calctiptiltdOTF(M5.*phase,p5,[3,3]);
                    seg5Prop = seg5(abs(seg5)>0);
                    seg5piston = mean(seg5Prop);
                    
                    
%                     R6 = sqrt((X-points{6}(2)).^2 + (Y-points{6}(1)).^2);
%                     M6 = R6<=mask_radius;
%                     seg6 = M6.*OPL;
%                     [tip6,tilt6] = calctiptiltdOTF(M6.*phase,p6,[3,3]);
%                     seg6Prop = seg6(abs(seg6)>0);
%                     seg6piston = mean(seg6Prop);
                    
                    
                    R7 = sqrt((X-points{7}(2)).^2 + (Y-points{7}(1)).^2);
                    M7 = R7<=mask_radius;
                    seg7 = M7.*OPL;
                    [tip7,tilt7] = calctiptiltdOTF(M7.*phase,p7,[3,3]);
                    seg7Prop = seg7(abs(seg7)>0);
                    seg7piston = mean(seg7Prop);
                    
                    
                    R8 = sqrt((X-points{8}(2)).^2 + (Y-points{8}(1)).^2);
                    M8 = R8<=mask_radius;
                    seg8 = M8.*OPL;
                    [tip8,tilt8] = calctiptiltdOTF(M8.*phase,p8,[3,3]);
                    seg8Prop = seg8(abs(seg8)>0);
                    seg8piston = mean(seg8Prop);
                    
                    
                    R9 = sqrt((X-points{9}(2)).^2 + (Y-points{9}(1)).^2);
                    M9 = R9<=mask_radius;
                    seg9 = M9.*OPL;
                    [tip9,tilt9] = calctiptiltdOTF(M9.*phase,p9,[3,3]);
                    seg9Prop = seg9(abs(seg9)>0);
                    seg9piston = mean(seg9Prop);
                    
                    
                    R10 = sqrt((X-points{10}(2)).^2 + (Y-points{10}(1)).^2);
                    M10 = R10<=mask_radius;
                    seg10 = M10.*OPL;
                    [tip10,tilt10] = calctiptiltdOTF(M10.*phase,p10,[3,3]);
                    seg10Prop = seg10(abs(seg10)>0);
                    seg10piston = mean(seg10Prop);
                    
                    
                    R11 = sqrt((X-points{11}(2)).^2 + (Y-points{11}(1)).^2);
                    M11 = R11<=mask_radius;
                    seg11 = M11.*OPL;
                    [tip11,tilt11] = calctiptiltdOTF(M11.*phase,p11,[3,3]);
                    seg11Prop = seg11(abs(seg11)>0);
                    seg11piston = mean(seg11Prop);
                    
                    
                    R12 = sqrt((X-points{12}(2)).^2 + (Y-points{12}(1)).^2);
                    M12 = R12<=mask_radius;
                    seg12 = M12.*OPL;
                    [tip12,tilt12] = calctiptiltdOTF(M12.*phase,p12,[3,3]);
                    seg12Prop = seg12(abs(seg12)>0);
                    seg12piston = mean(seg12Prop);
                    
                    
                    R13 = sqrt((X-points{13}(2)).^2 + (Y-points{13}(1)).^2);
                    M13 = R13<=mask_radius;
                    seg13 = M13.*OPL;
                    [tip13,tilt13] = calctiptiltdOTF(M13.*phase,p13,[3,3]);
                    seg13Prop = seg13(abs(seg13)>0);
                    seg13piston = mean(seg13Prop);
                    
                    
                    R14 = sqrt((X-points{14}(2)).^2 + (Y-points{14}(1)).^2);
                    M14 = R14<=mask_radius;
                    seg14 = M14.*OPL;
                    [tip14,tilt14] = calctiptiltdOTF(M14.*phase,p14,[3,3]);
                    seg14Prop = seg14(abs(seg14)>0);
                    seg14piston = mean(seg14Prop);
                    
                    
                    R15 = sqrt((X-points{15}(2)).^2 + (Y-points{15}(1)).^2);
                    M15 = R15<=mask_radius;
                    seg15 = M15.*OPL;
                    [tip15,tilt15] = calctiptiltdOTF(M15.*phase,p15,[3,3]);
                    seg15Prop = seg15(abs(seg15)>0);
                    seg15piston = mean(seg15Prop);
                    
                    
%                     R16 = sqrt((X-points{16}(2)).^2 + (Y-points{16}(1)).^2);
%                     M16 = R16<=mask_radius;
%                     seg16 = M16.*OPL;
%                     [tip16,tilt16] = calctiptiltdOTF(M16.*phase,p16,[3,3]);
%                     seg16Prop = seg16(abs(seg16)>0);
%                     seg16piston = mean(seg16Prop);
                    
                    
                    R17 = sqrt((X-points{17}(2)).^2 + (Y-points{17}(1)).^2);
                    M17 = R17<=mask_radius;
                    seg17 = M17.*OPL;
                    [tip17,tilt17] = calctiptiltdOTF(M17.*phase,p17,[3,3]);
                    seg17Prop = seg17(abs(seg17)>0);
                    seg17piston = mean(seg17Prop);
                    
                    
                    R18 = sqrt((X-points{18}(2)).^2 + (Y-points{18}(1)).^2);
                    M18 = R18<=mask_radius;
                    seg18 = M18.*OPL;
                    [tip18,tilt18] = calctiptiltdOTF(M18.*phase,p18,[3,3]);
                    seg18Prop = seg18(abs(seg18)>0);
                    seg18piston = mean(seg18Prop);
                    
                    R19 = sqrt((X-points{19}(2)).^2 + (Y-points{19}(1)).^2);
                    M19 = R19<=mask_radius;
                    seg19 = M19.*OPL;
                    [tip19,tilt19] = calctiptiltdOTF(M19.*phase,p19,[3,3]);
                    seg19Prop = seg19(abs(seg19)>0);
                    seg19piston = mean(seg19Prop);
                    
%                     R20 = sqrt((X-points{20}(2)).^2 + (Y-points{20}(1)).^2);
%                     M20 = R20<=mask_radius;
%                     seg20 = M20.*OPL;
%                     [tip20,tilt20] = calctiptiltdOTF(M20.*phase,p20,[3,3]);
%                     seg20Prop = seg20(abs(seg20)>0);
%                     seg20piston = mean(seg20Prop);
%                     
                    R21 = sqrt((X-points{21}(2)).^2 + (Y-points{21}(1)).^2);
                    M21 = R21<=mask_radius;
                    seg21 = M21.*OPL;
                    [tip21,tilt21] = calctiptiltdOTF(M21.*phase,p21,[3,3]);
                    seg21Prop = seg21(abs(seg21)>0);
                    seg21piston = mean(seg21Prop);
                    
                    R22 = sqrt((X-points{22}(2)).^2 + (Y-points{22}(1)).^2);
                    M22 = R22<=mask_radius;
                    seg22 = M22.*OPL;
                    [tip22,tilt22] = calctiptiltdOTF(M22.*phase,p22,[3,3]);
                    seg22Prop = seg22(abs(seg22)>0);
                    seg22piston = mean(seg22Prop);
                    
                    R23 = sqrt((X-points{23}(2)).^2 + (Y-points{23}(1)).^2);
                    M23 = R23<=mask_radius;
                    seg23 = M23.*OPL;
                    [tip23,tilt23] = calctiptiltdOTF(M23.*phase,p23,[3,3]);
                    seg23Prop = seg23(abs(seg23)>0);
                    seg23piston = mean(seg23Prop);
                    
                    R24 = sqrt((X-points{24}(2)).^2 + (Y-points{24}(1)).^2);
                    M24 = R24<=mask_radius;
                    seg24 = M24.*OPL;
                    [tip24,tilt24] = calctiptiltdOTF(M24.*phase,p24,[3,3]);
                    seg24Prop = seg24(abs(seg24)>0);
                    seg24piston = mean(seg24Prop);
                    
                    R25 = sqrt((X-points{25}(2)).^2 + (Y-points{25}(1)).^2);
                    M25 = R25<=mask_radius;
                    seg25 = M25.*OPL;
                    [tip25,tilt25] = calctiptiltdOTF(M25.*phase,p25,[3,3]);
                    seg25Prop = seg25(abs(seg25)>0);
                    seg25piston = mean(seg25Prop);
                    
                    R26 = sqrt((X-points{26}(2)).^2 + (Y-points{26}(1)).^2);
                    M26 = R26<=mask_radius;
                    seg26 = M26.*OPL;
                    [tip26,tilt26] = calctiptiltdOTF(M26.*phase,p26,[3,3]);
                    seg26Prop = seg26(abs(seg26)>0);
                    seg26piston = mean(seg26Prop);
                    
                    R27 = sqrt((X-points{27}(2)).^2 + (Y-points{27}(1)).^2);
                    M27 = R27<=mask_radius;
                    seg27 = M27.*OPL;
                    [tip27,tilt27] = calctiptiltdOTF(M27.*phase,p27,[3,3]);
                    seg27Prop = seg27(abs(seg27)>0);
                    seg27piston = mean(seg27Prop);
                    
                    R28 = sqrt((X-points{28}(2)).^2 + (Y-points{28}(1)).^2);
                    M28 = R28<=mask_radius;
                    seg28 = M28.*OPL;
                    [tip28,tilt28] = calctiptiltdOTF(M28.*phase,p28,[3,3]);
                    seg28Prop = seg28(abs(seg28)>0);
                    seg28piston = mean(seg28Prop);
%                     
%                     R29 = sqrt((X-points{29}(2)).^2 + (Y-points{29}(1)).^2);
%                     M29 = R29<=mask_radius;
%                     seg29 = M29.*OPL;
%                     [tip29,tilt29] = calctiptiltdOTF(M29.*phase,p29,[3,3]);
%                     seg29Prop = seg29(abs(seg29)>0);
%                     seg29piston = mean(seg29Prop);
%                     
%                     R30 = sqrt((X-points{30}(2)).^2 + (Y-points{30}(1)).^2);
%                     M30 = R30<=mask_radius;
%                     seg30 = M30.*OPL;
%                     [tip30,tilt30] = calctiptiltdOTF(M30.*phase,p30,[3,3]);
%                     seg30Prop = seg30(abs(seg30)>0);
%                     seg30piston = mean(seg30Prop);
%                     
%                     R31 = sqrt((X-points{31}(2)).^2 + (Y-points{31}(1)).^2);
%                     M31 = R31<=mask_radius;
%                     seg31 = M31.*OPL;
%                     [tip31,tilt31] = calctiptiltdOTF(M31.*phase,p31,[3,3]);
%                     seg31Prop = seg31(abs(seg31)>0);
%                     seg31piston = mean(seg31Prop);
%                     
%                     R32 = sqrt((X-points{32}(2)).^2 + (Y-points{32}(1)).^2);
%                     M32 = R32<=mask_radius;
%                     seg32 = M32.*OPL;
%                     [tip32,tilt32] = calctiptiltdOTF(M32.*phase,p32,[3,3]);
%                     seg32Prop = seg32(abs(seg32)>0);
%                     seg32piston = mean(seg32Prop);
%                     
%                     R33 = sqrt((X-points{33}(2)).^2 + (Y-points{33}(1)).^2);
%                     M33 = R33<=mask_radius;
%                     seg33 = M33.*OPL;
%                     [tip33,tilt33] = calctiptiltdOTF(M33.*phase,p33,[3,3]);
%                     seg33Prop = seg33(abs(seg33)>0);
%                     seg33piston = mean(seg33Prop);
%                     
%                     R34 = sqrt((X-points{34}(2)).^2 + (Y-points{34}(1)).^2);
%                     M34 = R34<=mask_radius;
%                     seg34 = M34.*OPL;
%                     [tip34,tilt34] = calctiptiltdOTF(M34.*phase,p34,[3,3]);
%                     seg34Prop = seg34(abs(seg34)>0);
%                     seg34piston = mean(seg34Prop);
%                     
%                     R35 = sqrt((X-points{35}(2)).^2 + (Y-points{35}(1)).^2);
%                     M35 = R35<=mask_radius;
%                     seg35 = M35.*OPL;
%                     [tip35,tilt35] = calctiptiltdOTF(M35.*phase,p35,[3,3]);
%                     seg35Prop = seg35(abs(seg35)>0);
%                     seg35piston = mean(seg35Prop);
                    
                    R36 = sqrt((X-points{36}(2)).^2 + (Y-points{36}(1)).^2);
                    M36 = R36<=mask_radius;
                    seg36 = M36.*OPL;
                    [tip36,tilt36] = calctiptiltdOTF(M36.*phase,p36,[3,3]);
                    seg36Prop = seg36(abs(seg36)>0);
                    seg36piston = mean(seg36Prop);
                    
                    R37 = sqrt((X-points{37}(2)).^2 + (Y-points{37}(1)).^2);
                    M37 = R37<=mask_radius;
                    seg37 = M37.*OPL;
                    [tip37,tilt37] = calctiptiltdOTF(M37.*phase,p37,[3,3]);
                    seg37Prop = seg37(abs(seg37)>0);
                    seg37piston = mean(seg37Prop);
                    
                    
                    
                    
                    PTTpos = zeros(37,3);
                    PTTpos(1,1) = seg1piston;
                    PTTpos(1,2) = tip1;
                    PTTpos(1,3) = tilt1;
                    PTTpos(2,1) = seg2piston;
                    PTTpos(2,2) = tip2;
                    PTTpos(2,3) = tilt2;
                    PTTpos(3,1) = seg3piston;
                    PTTpos(3,2) = tip3;
                    PTTpos(3,3) = tilt3;
                    PTTpos(4,1) = seg4piston;
                    PTTpos(4,2) = tip4;
                    PTTpos(4,3) = tilt4;
                    PTTpos(5,1) = seg5piston;
                    PTTpos(5,2) = tip5;
                    PTTpos(5,3) = tilt5;
%                     PTTpos(6,1) = seg6piston;
%                     PTTpos(6,2) = tip6;
%                     PTTpos(6,3) = tilt6;
                    PTTpos(7,1) = seg7piston;
                    PTTpos(7,2) = tip7;
                    PTTpos(7,3) = tilt7;
                    PTTpos(8,1) = seg8piston;
                    PTTpos(8,2) = tip8;
                    PTTpos(8,3) = tilt8;
                    PTTpos(9,1) = seg9piston;
                    PTTpos(9,2) = tip9;
                    PTTpos(9,3) = tilt9;
                    PTTpos(10,1) = seg10piston;
                    PTTpos(10,2) = tip10;
                    PTTpos(10,3) = tilt10;
                    PTTpos(11,1) = seg11piston;
                    PTTpos(11,2) = tip11;
                    PTTpos(11,3) = tilt11;
                    PTTpos(12,1) = seg12piston;
                    PTTpos(12,2) = tip12;
                    PTTpos(12,3) = tilt12;
                    PTTpos(13,1) = seg13piston;
                    PTTpos(13,2) = tip13;
                    PTTpos(13,3) = tilt13;
                    PTTpos(14,1) = seg14piston;
                    PTTpos(14,2) = tip14;
                    PTTpos(14,3) = tilt14;
                    PTTpos(15,1) = seg15piston;
                    PTTpos(15,2) = tip15;
                    PTTpos(15,3) = tilt15;
%                     PTTpos(16,1) = seg16piston;
%                     PTTpos(16,2) = tip16;
%                     PTTpos(16,3) = tilt16;
                    PTTpos(17,1) = seg17piston;
                    PTTpos(17,2) = tip17;
                    PTTpos(17,3) = tilt17;
                    PTTpos(18,1) = seg18piston;
                    PTTpos(18,2) = tip18;
                    PTTpos(18,3) = tilt18;
                    PTTpos(19,1) = seg19piston;
                    PTTpos(19,2) = tip19;
                    PTTpos(19,3) = tilt19;
%                     PTTpos(20,1) = seg20piston; %pokeseg
%                     PTTpos(20,2) = tip20;
%                     PTTpos(20,3) = tilt20;
                    PTTpos(21,1) = seg21piston;
                    PTTpos(21,2) = tip21;
                    PTTpos(21,3) = tilt21;
                    PTTpos(22,1) = seg22piston;
                    PTTpos(22,2) = tip22;
                    PTTpos(22,3) = tilt22;
                    PTTpos(23,1) = seg23piston;
                    PTTpos(23,2) = tip23;
                    PTTpos(23,3) = tilt23;
                    PTTpos(24,1) = seg24piston;
                    PTTpos(24,2) = tip24;
                    PTTpos(24,3) = tilt24;
                    PTTpos(25,1) = seg25piston;
                    PTTpos(25,2) = tip25;
                    PTTpos(25,3) = tilt25;
                    PTTpos(26,1) = seg26piston;
                    PTTpos(26,2) = tip26;
                    PTTpos(26,3) = tilt26;
                    PTTpos(27,1) = seg27piston;
                    PTTpos(27,2) = tip27;
                    PTTpos(27,3) = tilt27;
                    PTTpos(28,1) = seg28piston;
                    PTTpos(28,2) = tip28;
                    PTTpos(28,3) = tilt28;
%                     PTTpos(29,1) = seg29piston;
%                     PTTpos(29,2) = tip29;
%                     PTTpos(29,3) = tilt29;
%                     PTTpos(30,1) = seg30piston;
%                     PTTpos(30,2) = tip30;
%                     PTTpos(30,3) = tilt30;
%                     PTTpos(31,1) = seg31piston;
%                     PTTpos(31,2) = tip31;
%                     PTTpos(31,3) = tilt31;
%                     PTTpos(32,1) = seg32piston;
%                     PTTpos(32,2) = tip32;
%                     PTTpos(32,3) = tilt32;
%                     PTTpos(33,1) = seg33piston;
%                     PTTpos(33,2) = tip33;
%                     PTTpos(33,3) = tilt33;
%                     PTTpos(34,1) = seg34piston;
%                     PTTpos(34,2) = tip34;
%                     PTTpos(34,3) = tilt34;
%                     PTTpos(35,1) = seg35piston;
%                     PTTpos(35,2) = tip35;
%                     PTTpos(35,3) = tilt35;
                    PTTpos(36,1) = seg36piston;
                    PTTpos(36,2) = tip36;
                    PTTpos(36,3) = tilt36;
                    PTTpos(37,1) = seg37piston;
                    PTTpos(37,2) = tip37;
                    PTTpos(37,3) = tilt37;
                    
                    
                    
                    
                    

                    PTTpos
                    input('Press Enter to Send to DM');
                    
                    cd /home/lab/Desktop/Shared_Stuff
                    save('PTTpos','PTTpos');
                else
                    fprintf('Not Sending Data to DM\n');
                end
            end
    end
    
    if select == 2
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
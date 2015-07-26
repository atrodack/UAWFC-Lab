function CUBE = Automated_Testbed_run_cam_ver2(Run_Cam_Parameters)
%________________________________________________________________________
%Contact: K.Miller [millerk2@email.arizona.edu]               July 2015
%Contact: A.Rodack [atrodack@email.arizona.edu]
%Contact: J.Knight [jknight@optics.arizona.edu]
%
%Opens imaqtool and take takes some pictures. Input sets all controllable
%commands. Output are data structures for the Darks (if taken) and the PSFs
%_________________________________________________________________________
% INPUTS:
% Run_Cam_Parameters{}
% 1: Number of Triggers (0 based, want 10 triggers, this equals 9)
% 2: Number of Frames per Trigger
% 3: Shutter Speed of Camera
% 4: Take Darks: T/F
% 5: Camera Channel
% 8: Gain
% 9: Exposure
% METADATA Cells
% 6: Laser Type
% 7: Filter Type

%% Grab Information From Run_Cam_Parameters
ntriggers = Run_Cam_Parameters{1};
nframespertrigger = Run_Cam_Parameters{2};
sSpeed = Run_Cam_Parameters{3};
TakeDarks = Run_Cam_Parameters{4};
Channel = Run_Cam_Parameters{5};
Laser_Type = Run_Cam_Parameters{6};
Filter_Type = Run_Cam_Parameters{7};
Gain = Run_Cam_Parameters{8};
Exposure = Run_Cam_Parameters{9};

%% Initialize CUBE
CUBE = struct;
CUBE.PSFs = [];
CUBE.DARKS = [];
CUBE.BACKGROUND = [];
CUBE.Camera_Data = [];

nframes = nframespertrigger * (ntriggers+1);
fps = 30;
pausetime = nframes / fps;

%% Open imaqtool
imaqtool

%% Camera Acquistion and Settings
%Recognizes and labels FPM and Science Image and sources
scim = videoinput('dcam',Channel,'Y16_640x480'); %updated adaptorname to match matlab adaptorname from imaqhwinfo on 2/4/2015 by Justin

srcscim=getselectedsource(scim);

scim.FramesPerTrigger=nframespertrigger;
scim.TriggerRepeat = ntriggers;

srcscim.BusSpeed='S400';
srcscim.FrameRate='30'; 



pause(1);
preview(scim);
srcscim.Shutter = sSpeed;
srcscim.Gain = Gain;
srcscim.Exposure = Exposure;

%% Darks
if TakeDarks == true
   
%    input('Press Enter to Take Background');
%    closepreview;
%    start(scim);
%    BACKGROUND = getdata(scim,1);
%    BACKGROUND = int16(BACKGROUND / (2^4));
%    CUBE.BACKGROUND = BACKGROUND;
   
   fprintf('\n');
   input('Press Enter to take Darks');
   
       fprintf('Taking Frame %d\n',n);
       start(scim);
       pause(pausetime +2);
       
       dark = getdata(scim,nframes);
       dark = int16(dark / (2^4));
       CUBE.DARKS = dark;
   
    
   fprintf('Darks Complete\n');
   input('Press Enter to Take PSFs');
else
%     input('Press Enter to Close Preview');
%     closepreview;
    CUBE.DARKS = zeros(480,640,1,nframes);
    CUBE.BACKGROUND = 0;
   
end


%% PSFs
pause(5);

closepreview;


    start(scim);
    pause(pausetime +2);
    
    [scimage, ~, scimetadata] = getdata(scim,nframes);
    scimage = int16(scimage / (2^4));
    
    
    % Add some info into scimetadata(1) This will only be in the first
    % image data
    scimetadata(1).Shutter = sSpeed;
    if Channel == 1
        scimetadata(1).Channel = 'Science_Camera';
    elseif Channel == 2
        scimetadata(1).Channel = 'FPM_Camera';
    elseif Channel == 3
        scimetadata(1).Channel = 'Lyot_Stop_Camera';
    end
    scimetadata(1).Laser_Type = Laser_Type;
    scimetadata(1).Filter_Type = Filter_Type;
    
    CUBE.PSFs = scimage;


CUBE.Camera_Data = scimetadata;
    
    
















































































end

function CUBE = Testbed_run_cam(Run_Cam_Parameters)
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
% 1: Number of Images to Take
% 2: Shutter Speed of Camera
% 3: Take Darks: T/F
% 4: Camera Channel
% METADATA Cells
% 5: Laser Type
% 6: Filter Type

%% Grab Information From Run_Cam_Parameters
nframes = Run_Cam_Parameters{1};
sSpeed = Run_Cam_Parameters{2};
TakeDarks = Run_Cam_Parameters{3};
SciCamChannel = Run_Cam_Parameters{4};
Laser_Type = Run_Cam_Parameters{5};
Filter_Type = Run_Cam_Parameters{6};

%% Initialize CUBE
CUBE = struct;
CUBE.PSFs = [];
CUBE.DARKS = [];
CUBE.BACKGROUND = [];
CUBE.Camera_Data = [];
CUBE.Laser_Type = Laser_Type;
CUBE.Filter_Type = Filter_Type;


%% Open imaqtool
imaqtool

%% Camera Acquistion and Settings
%Recognizes and labels FPM and Science Image and sources
scim = videoinput('dcam',SciCamChannel,'Y16_640x480'); %updated adaptorname to match matlab adaptorname from imaqhwinfo on 2/4/2015 by Justin

srcscim=getselectedsource(scim);

%Sets camera settings
%src.GammaMode = 'manual';
%src.Gamma = 0;
scim.FramesPerTrigger=1;
srcscim.BusSpeed='S400';
srcscim.FrameRate='30'; 
srcscim.Shutter = sSpeed;



%% Darks
if TakeDarks == true
   
   input('Press Enter to Take Background');
   start(scim);
   BACKGROUND = getdata(scim,1);
   BACKGROUND = int16(BACKGROUND / (2^4));
   CUBE.BACKGROUND = BACKGROUND;
   
   fprintf('\n');
   input('Press Enter to take Darks');
   for n = 1:nframes
       fprintf('Taking Frame %d\n',n);
       start(scim);
       dark = getdata(scim,1);
       dark = int16(dark / (2^4));
       CUBE.DARKS(:,:,n) = dark;
   end
    
   fprintf('Darks Complete\n');
        
end


%% PSFs

input('Press Enter to Take PSFs');
metadata = cell(nframes,1);
PSFs = 0;
for n = 1:nframes
    fprintf('Taking Frame %d\n',n);
    start(scim);
    
    [scimage, scitime, scimetadata] = getdata(scim,1);
    scimage = int16(scimage / (2^4));
    metadata{n} = (scimetadata);
    CUBE.PSFs(:,:,n) = scimage;
end


CUBE.Camera_Data = metadata;
    
    
















































































end

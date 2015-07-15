

clear all;
clc;

%% Setup Folder Locations

directory1 = '/home/alex/Desktop/Data/NECO/';
var1 = 'dOTF_experiment_04_04_2013';
var2 = 'dOTF_experiment_04_11_2013';

directorypath1 = sprintf('%s/%s',directory1,var1);
directorypath2 = sprintf('%s/%s',directory1,var2);

foldername1 = cell(4,1);
foldername1{1} = 'Thu_Apr_04_Grab_Results_TiltAmp_.5_Seg_18';
foldername1{2} = 'Thu_Apr_04_Grab_Results_TiltAmp_.05_Seg_18';
foldername1{3} = 'Thu_Apr_04_Grab_Results_TiltAmp_.25_Seg_18';
foldername1{4} = 'Thu_Apr_04_Grab_Results_TiltAmp_1_Seg_18';

foldername2 = cell(4,1);
foldername2{1} = '04_11_2013_Beacon_Defocused_31mm_TiltAmp_1_Seg20';
foldername2{2} = '04_11_2013_Beacon_Focused_TiltAmp_1_Seg20';
foldername2{3} = '04_11_2013_Laser_Defocused_31mm_TiltAmp_1_Seg20';
foldername2{4} = '04_11_2013_Laser_Focused_TiltAmp_1_Seg20';


%% Choose Folder to Go To
chosen_path = directorypath2;
chosen_folder = foldername2{4};
% chosen_folder = foldername2{4};


current_folder = pwd;
cd(sprintf('%s/%s',chosen_path,chosen_folder));

%% Experimental Setup
lambda = AOField.HeNe_Laser;
pixelshift = [1,1];



%% Do Stuff with the Data
if chosen_path == directorypath1
    overlapsegs = [18];
    switch chosen_folder
        case foldername1{1}
            load('dOTFs.mat');
            load('VALS.mat');
            
            
            
        case foldername1{2}
            load('dOTFs.mat');
            load('VALS.mat');
            
            
            
        case foldername1{3}
            load('dOTFs.mat')
            
            
            
        case foldername1{4}
            load('dOTFs.mat');
            load('DISPLACEMENTS.mat');
            
            framespersec = 30;
            cubesize = size(dOCUBE);
            numframes = cubesize(3);
            
            movietime = numframes / framespersec;
            
%             for n = 1:numframes
%                 plotComplex(dOCUBE(:,:,n),3);
%                 bigtitle(sprintf('Time = %0.3f',n / framespersec),15);
%                 drawnow;
%                 pause(0.15);
%             end
            
            clf;
            
%             for n = 1:numframes
%                 plot(ALL_DISPLACEMENTS(:,1,n),ALL_DISPLACEMENTS(:,2,n),'r*')
%                 bigtitle(sprintf('SH Centroids at time = %0.3f',n / framespersec),15);
%                 drawnow;
%                 pause(0.15);
%             end
           
            closelooptime = 0.4; %frame number / framespersec
            closeloopframe = closelooptime * framespersec;
            
            wfsbias = 0;
            for n = closeloopframe:numframes
                wfsbias = wfsbias + ALL_DISPLACEMENTS(:,:,n);
            end
            wfsbias = wfsbias / (n - closeloopframe);
            
            clf;
            
%             for n = 1:numframes
%                 quiver(wfsbias(:,1),wfsbias(:,2),ALL_DISPLACEMENTS(:,1,n),ALL_DISPLACEMENTS(:,2,n))
%                 bigtitle(sprintf('Time = %0.3f',n / framespersec),15);
%                 drawnow;
%                 pause(0.15);
%             end
            
        otherwise
    end
    
    
elseif chosen_path == directorypath2
    overlapsegs = [20];
    switch chosen_folder
        case foldername2{1}
            load('dOTF_.mat'); %average dOTF from data set
            load('dOTFs_AminusB.mat');
            
            
            
            
        case foldername2{2}
            load('dOTF_.mat'); %average dOTF from data set
            load('dOTFs_AminusB.mat');
            
            
            
            
            
        case foldername2{3}
            load('dOTF_.mat'); %average dOTF from data set
            load('dOTFs_AminusB.mat');
            
            
            
            
            
            
        case foldername2{4}
            load('dOTF_.mat'); %average dOTF from data set
            load('dOTFs_AminusB.mat');
            
            
            
            
            
        otherwise
            
    end
end

% % load('dOTF_.mat'); %average dOTF from data set
% load('dOTFs.mat');
% % load('VALS.mat');
% load('DISPLACEMENTS.mat');


% cd(sprintf('%s/%s',directorypath2,foldername2{1}));
cd(current_folder)




    


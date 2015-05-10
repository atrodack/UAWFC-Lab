function [ read_in_images ] = BatchRead( Num_Folders, Num_Files_per_Folder, show_images, varargin)
%[ read_in_images ] = BatchRead( Num_Folders, Num_Files_per_Folder, varargin )
% Num_Folders is the number of Batch Folders to read in
% Num_Files_per_Folder is the number of files in the Batch
% varargin{1} through varargin{Num_Folders}: is path to folder
% varargin{Num_Folders+1} through varargin{Num_Folders * 2} is filename indicator

current_directory = cd;
varargout = cell(Num_Folders,1);

for n = 1:Num_Folders
    
    cd(char(varargin{1}(n)));
    
    for m = 1:Num_Files_per_Folder
        filename = sprintf('%s%d.fits',char(varargin{1}(Num_Folders+n)),m);
        pic = fitsread(filename);
        
        if show_images == true
            clf;
            figure(1)
            imagesc(pic);
            title(sprintf('Reading in Picture %d',m));
            sqar;
            colormap(gray);
            drawnow;
            pause(1);
        end
        
        img(:,:,m) = double(pic);
    end
    
    varargout{n} = img;
end
read_in_images = varargout;
cd(current_directory);
end


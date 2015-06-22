function [CUBE,FILES] = cubifyImages(REGEXP,DARK,PLANE)

% function [CUBE,FILES] = cubifyImages(REGEXP,[DARK],[PLANE])
% 
% 20120416 jlc: Added better support for varying image sizes.
% 20120504 jlc: incorporated dark subtraction.

FILES = dir(REGEXP);

if(nargin<2)
    DARK = 0;
end

if(nargin<3)
    PLANE = 1;
end

DARK = int16(DARK);

% SIZES = zeros([length(FILES), 2]);
% for n=1:length(FILES)
%     finfo = imfinfo(FILES(n).name);
%     SIZES(n,1) = finfo.Width;
%     SIZES(n,2) = finfo.Height;
% end

INFO1 = imfinfo(FILES(1).name);
IMG = (imread(FILES(1).name));

% SIZE = 2.^(ceil(log2(max(SIZES))));
% CUBE = zeros([INFO1.Width, INFO1.Height,length(FILES)],'int16'); % WARNING: possible loss of top bit?
CUBE = zeros([INFO1.Height, INFO1.Width,length(FILES)],class(IMG)); % WARNING: possible loss of top bit?

for n=1:length(FILES)
    modprint(n,25);
    %INFO = imfinfo(FILES(n).name);
    %IMG = int16(imread(FILES(n).name));
    IMG = (imread(FILES(n).name));
    IMG = IMG(:,:,PLANE);
    %CH1 = 1:size(IMG,1);
    %CH2 = 1:size(IMG,2);
    %MIN = (IMG(1,1)+IMG(1,end)+IMG(end,1)+IMG(end,end))/4;
    %CUBE(:,:,n) = MIN;
    %CUBE(CH1,CH2,n) = IMG(:,:,1);
    %CUBE(CH1,CH2,n) = IMG;
    CUBE(:,:,n) = IMG;
end
fprintf('\n');

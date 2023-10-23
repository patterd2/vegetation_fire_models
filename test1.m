nImage = 150;                                                        % L1
fps = 3.0;                                                           % L2
addToInFolder = 'Address\to\Input\Images\Folder';                    % L3
addToOutFolder = 'Address\to\Output\Video\Folder';                   % L4

oVideo = VideoWriter(fullfile(addToOutFolder, 'myVideo.avi'));       % L5
oVideo.FrameRate = fps;                                              % L6
open(oVideo)                                                         % L7
for i = 1:nImage                                                     % L8
    fname = ['image' num2str(i, '%.2d') '.png'];                     % L9
    curImage = imread(fullfile(addToInFolder,fname));                % L10
    writeVideo(oVideo, curImage);                                    % L11
end                                                                  % L12
close(oVideo)                                                        % L13


Test = zeros(2,2);
Test(1,1) = 0; 
Test(2,1) = 1;
Test(1,2) = 2;
Test(2,2) = 0;
map = [0.4660 0.6740 0.1880
    0.0039 0.1953 0.1250
    0.8500 0.3250 0.0980
    0.5 0.5 0.5];
colormap(map)

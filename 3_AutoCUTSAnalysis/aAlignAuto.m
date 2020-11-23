% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
VISUALIZE = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read files %%%%%%%%%%%%%%%%%%%%%%%%%%%
name =  'Example';
folder1 ='AutuCUTS_Pipeline';
folder2 =[name,'_1_order'];
fileName = '1'; % Rember to rename all files to 1, so they are chronological ordered
ImgType = 'tif';
buffer = 0; % Define which image number to start from
factor = 4; % factor reduce the images
markBackground = 236; % Define background value

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create save image folder %%%%%%%%%%
[fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1);

% Get path to read files
[srcFiles,path]=pathFilesFnc(folder1,folder2,ImgType);
% Define first image
[fixed, markBackground] =firstImgInfoFnc(path,fileName,folderNameSave,fileNameSave,markBackground);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read Images and Align%%%%%%%%%%%%%%%%%%%%%%
if VISUALIZE
    h = gca;
end

for i =2:size(srcFiles,1)+buffer
    %%
    disp(['Register image number ',num2str(i),' out of ',num2str(length(srcFiles)+buffer)])
    K =imread([path,fileName,' (',num2str((i)),').tif']);
    moving = medfilt2(rgb2gray(K)); % Apply Median filter
    [moving_reg] =alignFnc(fixed,moving,factor,markBackground);
    if VISUALIZE
        imshow(uint8(moving_reg),'Parent',h); title(i);drawnow;
    end
    %% Save files
    imwrite(moving_reg, [folderNameSave,'\',fileNameSave,'_0',num2str(i),'.png'],'png');
    fixed = moving_reg;
end
%% %%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1)
% Create the save folder of images
% Sintax:
%     [fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1)
% Inputs:
%     name,        name of subject
%     folder1,     First folder name

% Outputs:
%     fileNameSave,         Output image name
%     folderNameSave,       Directory for saved folder

s = what(folder1);
savePath=s.path;
folderSave =[name,'_2_align'];
folder = [savePath,'\',folderSave];
fileNameSave ='Img_align';
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end

end

function [fixed, markBackground] =firstImgInfoFnc(path,fileName,folderNameSave,fileNameSave,markBackground)
% Read first image and define background, so the black color of alignment
% will get removed
% Sintax:
%     [fixed, markBackground] =firstImgInfoFnc(path,fileName,folderNameSave,fileNameSave,markBackground)
% Inputs:
%     path,               path directory of where images should be read from
%     fileName,           Read files with only this name in front
%     folderNameSave,     path directory of where to save aligned images
%     fileNameSave,       Output image name
%     markBackground      Intensity value of grayscale images that is
%                         regarded as background

% Outputs:
%     fixed,              Fixed image used for registration
%     markBackground      Intensity value of grayscale images that is
%                         regarded as background

try
    original = rgb2gray(imread([path,fileName,' (1).tif']));
    imwrite(medfilt2(original), [folderNameSave,'\',fileNameSave,'_01','.png'],'png');
    fixed = medfilt2(original);% Apply Median filter
    if isempty(markBackground)==1
        imshow(fixed); title('Press on background outside ROI');
        [x1,y1]=ginput(1);
        markBackground =fixed(round(x1),round(y1));
        close all
    end
catch
    disp('Rename files in folder to 1, so they are chronological ordered')
end
end

function [srcFiles,path]=pathFilesFnc(folder1,folder2,ImgType)
% Get path to read files
% Sintax:
%     [srcFiles,path]=pathFilesFnc(folder1,folder2,ImgType)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     ImgType,     Image type

% Outputs:
%     srcFiles,     Structure with file informations
%     path,         path directory of where images should be read from

srcFiles = dir([fullfile(folder1,'/',folder2),'/*.',ImgType]);
path=srcFiles.folder;
path=[path '\'];

if isempty(srcFiles)==1
    folder = uigetdir;
    srcFiles=dir([fullfile(folder),'/*.',ImgType]);
    path=srcFiles.folder;
    path=[path '\'];
end
end


function [moving_reg] =alignFnc(fixed,moving,factor,markBackground)
% Align images with rigid registration
% Sintax:
%     [moving_reg] =alignFnc(fixed,moving,factor,markBackground)
% Inputs:
%     fixed,              Fixed image used for registration
%     moving,             Moved image that need to be registered 
%     factor,             How much the images are scaled down to increase
%                         the speed
%     markBackground      Intensity value of grayscale images that is
%                         regarded as background

% Outputs:
%     moving_reg          Registered image of the moving image
%  

disp('<strong>Intensity based alignment</strong>');
lastwarn('');
% Increase contrast and blur image of fixed and moving
% image
FIXED = imgaussfilt(adapthisteq(imresize(fixed, 1/factor)),1.3);
MOVING = imgaussfilt(adapthisteq(imresize(moving, 1/factor)),1.3);

% Default spatial referencing objects
fixedRefObj = imref2d(size(FIXED));
movingRefObj = imref2d(size(MOVING));

% Intensity-based registration
[optimizer, metric] = imregconfig('multimodal');

% Align centers
fixedCenterXWorld = mean(fixedRefObj.XWorldLimits);
fixedCenterYWorld = mean(fixedRefObj.YWorldLimits);
movingCenterXWorld = mean(movingRefObj.XWorldLimits);
movingCenterYWorld = mean(movingRefObj.YWorldLimits);
translationX = fixedCenterXWorld - movingCenterXWorld;
translationY = fixedCenterYWorld - movingCenterYWorld;

% Coarse alignment
initTform = affine2d();
initTform.T(3,1:2) = [translationX, translationY];

tform = imregtform(MOVING,movingRefObj,FIXED,fixedRefObj,'rigid',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);

% Apply transformation
tform.T(3,1:2) =tform.T(3,1:2)*factor;
outputView = imref2d(size(fixed));
movingRegistered  = imwarp(moving,tform,'OutputView',outputView);
movingRegistered(movingRegistered<4)=markBackground;
%% Second alignment
disp('Scond alignment');
%%%%%%%%%%%%%%%%%%%%Downscale registered image
MOVING = imgaussfilt(adapthisteq(imresize(movingRegistered, 1/factor)),1.3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Version 2
% Default spatial referencing objects
fixedRefObj = imref2d(size(FIXED));
movingRefObj = imref2d(size(MOVING));

% Intensity-based registration
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
metric.UseAllPixels = true;
optimizer.GrowthFactor = 1.050000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 1.35000e-03;
optimizer.MaximumIterations = 300;


% Align centers
fixedCenterXWorld = mean(fixedRefObj.XWorldLimits);
fixedCenterYWorld = mean(fixedRefObj.YWorldLimits);
movingCenterXWorld = mean(movingRefObj.XWorldLimits);
movingCenterYWorld = mean(movingRefObj.YWorldLimits);
translationX = fixedCenterXWorld - movingCenterXWorld;
translationY = fixedCenterYWorld - movingCenterYWorld;
% Coarse alignment
initTform = affine2d();
initTform.T(3,1:2) = [translationX, translationY];

% Apply Gaussian blur
fixedInit = imgaussfilt(FIXED,1.000000);
movingInit = imgaussfilt(MOVING,1.000000);

% Normalize images
movingInit = mat2gray(movingInit);
fixedInit = mat2gray(fixedInit);

% Apply transformation
tform2 = imregtform(movingInit,movingRefObj,fixedInit,fixedRefObj,'rigid',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
tform2.T(3,1:2) =tform2.T(3,1:2)*factor;

moving_reg = imwarp(movingRegistered,tform2,'OutputView',imref2d(size(fixed)));
moving_reg(moving_reg<4)=markBackground;
end
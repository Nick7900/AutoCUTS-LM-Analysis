% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
folder1 = 'Segmentation';
folder2 ='Example';  % Direct to folder where the images are
ImgType = 'tif'; % Image type
% Read files in directory
[folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType);
%%
for i = 1:length(srcFiles)
    try
        % Filename
        [~,name] = fileparts([srcFiles(i).folder,'\',srcFiles(i).name]);
        disp(['Process image number ',num2str(i),' out of ',num2str(length(srcFiles))])
        rgbImage = imread(strcat(folder1,'/',folder2,'/',srcFiles(i).name));          % Read image
        % Creating the processbar
        f = waitbar(0,'Loading Data','Name','Extract sections');
        pause(.1)
        % Run function
        [Data,f] =imgSegmentationFnc(rgbImage,f);
        waitbar(1,f,'Finishing','Name','Extract sections');
        pause(0.1)
        close(f)
        % Rename output to insure they are in chronological order
        [name] =renameOutputFnc(name,srcFiles,i);
        [folderSave]=folderGenerateFnc(folder1,folder2,1);
        % Save images
        for j = 1:length(Data)
            imwrite(Data(j).Images, [folder1,'/',folderSave,'/',name,'_0',num2str(Data(j).Index),'.tif'],'tif');
        end 
    catch ME
        %%
        disp(['mistake in picture: ',name,' #',num2str(i)]);
        [folderSave]=folderGenerateFnc(folder1,folder2,0);
        imwrite(rgbImage, [folder1,'/',folderSave,'/',name,'.',ImgType],ImgType);
        waitbar(1,f,'Finishing','Name','Extract sections');
        pause(1)
        close(f)
    end 
end
disp('Finished')

%% Function to get the Images
function [folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType)
% Get path to read files
% Sintax:
%     [folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     ImgType,     Image type (eg. TIF, png, JPEG)

% Outputs:
%     srcFiles,     Structure with file informations
%     folder2,      Second folder name
try
    srcFiles = dir([fullfile(folder1,'/',folder2),'/*.',ImgType]);

catch
    % Open diolog box if the folder is placed wrong
    if isempty(srcFiles)==1
        folder2 = uigetdir;
        srcFiles=dir([fullfile(folder2),'/*.',ImgType]);

    end
end
end

function  [ImgInfo,f] =imgSegmentationFnc(rgbImage,f)
% Segmentation of sections, ordered from left to right
% Sintax:
%     [ImgInfo,f] =imgSegmentationFnc(rgbImage,f)
% Inputs:
%     rgbImage,     input image
%     f,            process bar
%
% Outputs:
%     ImgInfo,   segmented images
%     f,         process bar

% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[~, ~, numberOfColorChannels] = size(rgbImage);
if numberOfColorChannels > 1
    % It's not really gray scale like we expected - it's color.
    % Convert to gray image
    grayImage = uint8(rgb2gray(rgbImage)); % Take blue channel.
else
    grayImage = uint8(rgbImage); % It's already gray scale.
end

%% Converting the image
% Convert values below 85 into background
grayImage(grayImage<85)=255;
% Apply gaussian blur
Iblur = imgaussfilt(grayImage,0.8);
% Apply entrophy filter
E = entropyfilt(Iblur);
Eim = rescale(E);

% Convert to binary image
BW1 = imbinarize(Eim, .7);
% Fill eventually holes of detected mask
binaryHoles = imfill(BW1,'holes');

filtSmall = 100000; % Define minimum pixel required, depends on resolution
% Remove small objects
binaryRemove = bwareaopen(binaryHoles, filtSmall);
%%
%%%%% Label sections
% Lavel each section to an integer-valued image where all pixels
% in the sections have values of 1, or 2, or 3, or ... etc.
labeledImage = bwlabel(binaryRemove, 8);     % Label each blob so we can make measurements of it
labeledImage =int8(labeledImage);

% Bar status
try
    waitbar(.33,f,'Detect sections','Name','Extract sections');
    pause(.1)
catch 
    f = waitbar(0.33,'Detect sections','Name','Extract sections');
    pause(.5)
end
% Get blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labeledImage, grayImage,'Area','BoundingBox');

%% More than one section then we need to process

% Area new scanner
allSectionAreas = [blobMeasurements.Area];
sectionArea=max([blobMeasurements.Area]);

areaSmall = sectionArea*0.30; % Size of objects that will be removed
remove_section = allSectionAreas < areaSmall;% find small objects.
logigalIdx1 = ones(1,length(blobMeasurements));
smallIdxRemove =find(remove_section & logigalIdx1);
bigImageIdx = find(logigalIdx1);
areaRemove =sort(smallIdxRemove);
smallIdxRemove1 =flip(areaRemove);

% Remove all small objects from blobmeasurement
for k = 1 : length(smallIdxRemove1)
    blobMeasurements(smallIdxRemove1(k))=[];
    allSectionAreas(smallIdxRemove1(k))=[];
end

%% See if there is more than more section in the whole image
if size(blobMeasurements,1)>1==1
    % Flip the order, so it runs from left to right
    blobMeasurements = flip(blobMeasurements);
    
    % Now we need to fine the number of sections for the big images
    % Extract the big images
    % Save individual images
    waitbar(.67,f,'Processing data','Name','Extract sections');
    pause(0.5)
else
    bigImageIdx = uint8(1);
    waitbar(.67,f,'Processing data','Name','Extract sections');
    pause(0.5)
end
%% Now use the BoundingBox as a mask on the original image.

if ~isempty(bigImageIdx)==1
    DataSingle(size(bigImageIdx,2)) = struct('Images',[],'Index',[]);
    
    for k = 1 : size(bigImageIdx,2)           % Loop through all blobs.
        % Find the bounding box of each section.
        thisBlobsBoundingBox = blobMeasurements(k).BoundingBox;  % Get list of pixels in current blob.
        
        % Extract out this section into it's own image.
        subImage = imcrop(rgbImage, thisBlobsBoundingBox);
        DataSingle(k).Images=subImage;
        DataSingle(k).Index=bigImageIdx(k);
    end
else
    DataSingle =[];
end

ImgInfo = DataSingle;
% %% Sort the struct
Afields = fieldnames(ImgInfo);
Acell = struct2cell(ImgInfo);
sz = size(Acell);
% Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)
% Make each field a column
Acell = Acell';                         % (MxN)xP
Acell = sortrows(Acell, 2);% Sort by first field "Idexes"
Acell = reshape(Acell', sz);% Put back into original cell array format
ImgInfo = cell2struct(Acell, Afields, 1);% Convert to Struct

%% Visualize the sections
for k = 1:length(ImgInfo)
    textFontSize = 13;
    if length(ImgInfo) <= 4
        %         Display the image with informative caption.
        subplot(2, 2, k);
        imshow(ImgInfo(k).Images);
        caption = sprintf('Section #%d', ...
            k);
        title(caption, 'FontSize', textFontSize);
    elseif length(ImgInfo) <= 6
        %         Display the image with informative caption.
        subplot(2, 3, k);
        imshow(ImgInfo(k).Images);
        caption = sprintf('Section #%d', ...
            k);
        title(caption, 'FontSize', textFontSize);
    else
        %         Display the image with informative caption.
        subplot(3, 3, k);
        imshow(ImgInfo(k).Images);
        caption = sprintf('Section #%d',k);
        title(caption, 'FontSize', textFontSize);
    end

end

end

function name =renameOutputFnc(name,srcFiles,i)
% Change part of name if it contains '_' to keep them in chronological
% order
% Sintax:
%     name =renameOutputFnc(name,srcFiles,i)
% Inputs:
%     name,         name of output file
%     srcFiles,     Structure with file informations
%     i,            Image number from loop

% Outputs:
%     name,         name of output file adjusted if necessary

% Find the underscoreLocation of the name
underscoreLocation = strfind(srcFiles(i).name , '_00');
if isempty(underscoreLocation)==0
    name(underscoreLocation) ='-';
    name([underscoreLocation+1 underscoreLocation+2]) =[];
end
underscoreLocation2 = strfind(srcFiles(i).name , '_01');
if isempty(underscoreLocation2)==0 && length(name)>9==1
    name(underscoreLocation2(1)) ='-';
    name(underscoreLocation2(1)+1) =[];
end
underscoreLocation2 = strfind(srcFiles(i).name , '_02');
if isempty(underscoreLocation2)==0 && length(name)>9==1
    name(underscoreLocation2(1)) ='-';
    name(underscoreLocation2(1)+1) =[];
end
end

function [folderSave]=folderGenerateFnc(folder1,folder2,flag)
% Create the save folder of images
% Sintax:
%     [folderSave]=folderGenerateFnc(folder1,folder2)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     flag,        Notice if the image will be saved in the Save or
%                  Fail folder (1=Save, 0=Fail)
% Outputs:
%     folderSave,         Folder name for saved images

if flag==1
    s = what(folder1);
    savePath=s.path;
    folderSave =[folder2,'_save'];
    folderName = [savePath,'\',folderSave];
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
else
    s = what(folder1);
    savePath=s.path;
    folderSave =[folder2,'_fail'];
    folderName = [savePath,'\',folderSave];
    if ~exist(folderName, 'dir')
        mkdir(folderName)
    end
    
end
end

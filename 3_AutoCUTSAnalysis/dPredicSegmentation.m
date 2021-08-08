% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% Analysing pyramidal cells
VISUALIZE = true;
SMOOTH =true; % smooth the 3D-reconstruction image
CROP = false;
PROFILE =false; % Analyse every single profiles
SAVEFILES = true;
SAVEIMAGES =true;
INTERPZ=false;
res = 272   ; % output resolution (nm) in xy-direction
xy = 272; % image resolution (nm) in xy-direction
z = 800; % z-axis resolution (nm)
buffer =0;
side =0; % Pial surface direction: left side 0, right side 1
%%%%%%%%%%%%%%%%%% Path for files
name = 'Example';
folder1 ='AutuCUTS_Pipeline';
folder2 =[name,'_4_Prediction'];
folder3 =[name,'_UNetDense_stack_layer3'];
nameID = 'subject1';
fileName = [nameID,'layer3_'];
ImgType = 'png';

%% Start program
% Define folders
[folder1,folder2,folder3]=fileInforFnc(folder1,folder2,folder3,ImgType);
% Run the program
runProgram(folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,name,side,VISUALIZE,SMOOTH,CROP,PROFILE,INTERPZ,SAVEFILES,SAVEIMAGES);

%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
function [folder1,folder2,folder3,srcFiles,path]=fileInforFnc(folder1,folder2,folder3,ImgType)
%%
% Get path to read files
% Sintax:
%     [folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     folder3,     Third folder name
%     ImgType,     Image type (eg. TIF, png, JPEG)

% Outputs:
%     folder1,      First folder name
%     folder2,      Second folder name
%     folder3,      Third folder name
%     srcFiles,     Structure with file informations
%     path,         Filepath directory
try
    srcFiles = dir([fullfile(folder1,'/',folder2,'/',folder3),'/*.',ImgType]);
    path=srcFiles.folder;
    path=[path '\'];
    if isempty(srcFiles)==1
        path = uigetdir;
        srcFiles=dir([fullfile(path),'/*.',ImgType]);
        folder1 = path;
        folder2 = path;
        folder3 = path;
        dashLocation = strfind(srcFiles(1).folder , '\');
        folder1(dashLocation(end-1):end)=[]; % remove last two sub folder
        folder1(1:dashLocation(end-2))=[]; % remove begining to third last sub folder
        folder2(dashLocation(end):end)=[]; % remove last sub folder
        folder2(1:dashLocation(end-1))=[]; % remove begining to 2nd last sub folder
        folder3(1:(dashLocation(end)))=[]; % remove begining to last sub folder
        path=[path '\'];
    end
catch
    % Open diolog box if the folder is placed wrong
    if isempty(srcFiles)==1
        path = uigetdir;
        srcFiles=dir([fullfile(path),'/*.',ImgType]);
        folder1 = path;
        folder2 = path;
        folder3 = path;
        dashLocation = strfind(srcFiles(1).folder , '\');
        folder1(dashLocation(end-1):end)=[]; % remove last two sub folder
        folder1(1:dashLocation(end-2))=[]; % remove begining to third last sub folder
        folder2(dashLocation(end):end)=[]; % remove last sub folder
        folder2(1:dashLocation(end-1))=[]; % remove begining to 2nd last sub folder
        folder3(1:(dashLocation(end)))=[]; % remove begining to last sub folder
        path=[path '\'];
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
function  runProgram(folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,name,side,VISUALIZE,SMOOTH,CROP,PROFILE,INTERPZ,SAVEFILES,SAVEIMAGES)
% Read one image of the time from folder and crop it
% Note - The computer was out of memory but reading and cropping the
% whole stack, which is the reason for this approach.
% Sintax:
%     [x] = readImgSingleLoopFnc(folder1,folder2,folder3,fileName,ImgType,cRect,i)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     folder3,     Second folder name
%     fileName,    Name of image that we read from
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     res,         Output voxel resolution
%     xy,          xy pixel resolution of images
%     z,           z resolution (distance between each image)
%     buffer,      starting point from analysis the images (eg. buffer 20
%                  => analysise from image 20)
%     name,        name of subject/folder that will be analysed
%     side,        Pial surface direction: left side 0, right side 1
%     VISUALIZE,   Visualize output =1,
%     SMOOTH,      Smooth 3D-reconstruction images
%     CROP,        Crop imagestack (True or false)
%     INTERPZ,     Define if the imagestack should have the same resolution
%                  as x and y axis
%     SAVEFILES,   Save data
%     SAVEIMAGES,  Save filtered images

%%%%%%%%%%%%%%%%% Start program
radius = 0; % filtering noise objects away.
[BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,~] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius,INTERPZ);
%%%%%%%%%%%% 3D analysis
% Define window of image stack
window = [size(BW,2).*res/1000, size(BW,1).*res/1000,size(BW,3).*z/1000];  %[xLeft, yTop, length, width, hieght].
[blobfilt,CCfilt] =blobFnc(BW,z,res,side);
%%%%%%%%%%%% 2D analysis
[CC2D,bloblFilit2D] =Analysis2D2(BW,blobfilt,CCfilt,res,side);

%%%%%%%%%%%%%%% Run Kmeans
disp('Kmeans 3D')
[kmeans3D]=kmeans3DFnc(window,blobfilt,bloblFilit2D,folder,res,z);
disp('Kmeans 2D')
[kmeans2D] =kmeans2DFnc(bloblFilit2D,kmeans3D);
if PROFILE
    %%%%%%%%%%%% Profile analysis
    [blobfiltProfile] =AnalysisAllProfiles2D(BW,blobfilt,CCfilt,res);
    [kmeans2DP] =kmeans2DFnc(blobfiltProfile,kmeans3D);
    saveData2DFnc(blobfiltProfile,kmeans2DP,folder)
end
%% Visualize and save files
if VISUALIZE == 0
    close all
    [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,kmeans3D,res,xy,z,folder,name,3,SAVEIMAGES);
    [pyramideNeuronImage2D,roundNeuronImage2D] =saveDataReconstructionFnc(bloblFilit2D,CC2D,BW,kmeans2D,res,xy,z,folder,name,2,SAVEIMAGES);
    if SAVEIMAGES
        %%%%%%%%%%%%%%%%%%% Save images of binary segmentation
        for i = 1:size(pyramideNeuronImage,3)
            disp(['Save Image number ',num2str(i)])
            imwrite(pyramideNeuronImage(:,:,i),[folderImageSave,'/',fileNameSave,num2str(i),'.png'],'png');
            imwrite(roundNeuronImage(:,:,i),[folderImageSave2,'/',fileNameSave2,num2str(i),'.png'],'png');
            imwrite(pyramideNeuronImage2D(:,:,i),[folderImageSave,'/',fileNameSave,num2str(i),'.png'],'png');
            imwrite(roundNeuronImage2D(:,:,i),[folderImageSave2,'/',fileNameSave2,num2str(i),'.png'],'png');
        end
        disp('All images have been saved')
    end
else
    %%
    disp('Visualize 3D-reconstruction of pyramidal cells')
    [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,bloblFilit2D,CC2D,BW,kmeans3D,res,z,folder,name,SAVEFILES,SAVEIMAGES,SMOOTH);
    if SAVEFILES
        %%%%%%%%%%%%%%%%%%% Save images of binary segmentation
        for i = 1:size(pyramideNeuronImage,3)
            disp(['Save Image number ',num2str(i)])
            imwrite(pyramideNeuronImage(:,:,i),[folderImageSave,'/',fileNameSave,num2str(i),'.png'],'png');
            imwrite(roundNeuronImage(:,:,i),[folderImageSave2,'/',fileNameSave2,num2str(i),'.png'],'png');
        end
        disp('All images have been saved')
    end
end
end

function [blobfilt,CCfilt] =blobFnc(BW,z,res,side)
%%
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobfilt,CCfilt] =blobFnc(BW,z,res)
% Inputs:
%     BW,           Binary image stack
%     z,            z pixel resolution of images
%     res,          Output voxel resolution
%     side,         Pial surface direction: left side 0, right side 1
% Outputs:
%     blobfilt,     Table with information of detected objects
%     CCfilt,       Filtered Connected components, returned as a structure

disp('Detect connections')
CC = bwconncomp(BW);
blobinfo = regionprops3(CC,'BoundingBox','Centroid','Volume','Orientation','VoxelList','PrincipalAxisLength','SurfaceArea');
[blobRemove, CCfilt]= blobFilt(blobinfo,CC,BW); % Remove cells that are only in 1 plane
window = [size(BW,2), size(BW,1),round(size(BW,3).*z/res)];  %[xLeft, yTop, length, width, hieght].
[blobfilt,CCfilt] = blobParameters(blobRemove,CCfilt,z,res,side,window);
end

function [BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,cRect] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius,INTERPZ)
% Get path to read files
% Sintax:
%     [BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,cRect] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius)
% Inputs:
%     CROP,             Crop stack of images (True or False)
%     folder1,          First folder name
%     folder2,          Second folder name
%     folder3,          Third folder name
%     fileName,         Name of image that we read from
%     ImgType,          Image type (eg. TIF, png, JPEG)
%     res,              Output voxel resolution
%     xy,               xy pixel resolution of images
%     z,                z resolution (distance between each image)
%     buffer,           starting point from analysis the images (eg. buffer 20
%                       => analysise from image 20)
%     radius,           Artefact size based on a defined radius of a perfect sphere.
%     interpZ,          Interpolate z-axis, so the voxel resolution is the same

% Outputs:
%     BW,               Binary image stack
%     folderImageSave,  Folder location for saved images of pyramidal cells
%     fileNameSave,     Name of saved images of pyramidal cells
%     folderImageSave2, Folder location for saved images of removed objects
%     folder,           Folder location for all results
%     cRect,            Defined imag area of cropped imagestack


%%%%%%%%%%%%%%%%%% Sitting variables

%%%%%%%%%%%% filter small cells away
% Filter noise from binary image stack with a defined radius
% Default radius =0;
filtVoxels =filtEstimateFnc(radius,res,z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create save image folder %%%%%%%%%%
%%
s = what(folder1);
savePath=[s.path,'\',folder2];
folderResults =[folder3,'_results'];
imgFilt =('ImageFilt');
folder = [savePath,'\',folderResults];
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end

folderImageSave = [savePath,'\',folderResults,'\',imgFilt];
fileNameSave =[fileName,'Filt'];
if ~exist(folderImageSave, 'dir')
    mkdir(folderImageSave);
end


imgFilt =('ImagRound');
folder = [savePath,'\',folderResults];
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end


folderImageSave2 = [savePath,'\',folderResults,'\',imgFilt];
fileNameSave2 =[fileName,'Round'];
if ~exist(folderImageSave2, 'dir')
    mkdir(folderImageSave2);
end


% Read Image
[I,Ic,cRect] = readImgFnc(folder1,folder2,folder3,fileName,ImgType,res,xy,buffer,CROP);
%%
% Filter and fill the holes of the images
if filtVoxels>0
    if ~isempty(Ic) == 1
        disp('Filter Crop image')
        [BW]= filtFillFnc(Ic, filtVoxels);
    else
        disp('Filter Original image')
        [BW]= filtFillFnc(I, filtVoxels);
    end
else
    
    if ~isempty(Ic) == 1
        BW=Ic;
    else
        BW=I;
    end
end

%%
clear I
if INTERPZ
    BW =interpImage(BW,xy,res,z);
end


end


function BW_stack =interpImage(BW,xy,res,z)

zStep =round(linspace(1,size(BW,3),round(size(BW,3)/40)));
xC =round(size(BW,2)*xy/res); % Convert image to desired resolution
yC =round(size(BW,1)*xy/res); % Convert 272 to 1000 nm pr pixel value

xL =linspace(1,size(BW,2),xC);
yL =linspace(1,size(BW,1),yC);

BW_stack = ([]);
for ii = 1:size(zStep,2)-1
    % dont change z
    zC =round(length(zStep(ii):zStep(ii+1))*z/res); % Convert image to desired resolution
    [X,Y,Z] = meshgrid(1:1:size(BW,2),1:1:size(BW,1),1:1:length(zStep(ii):zStep(ii+1))); 
    zL =linspace(1,length(zStep(ii):zStep(ii+1)),zC);
    [Xq,Yq,Zq] = meshgrid(xL,yL,zL);
    I_interp = interp3(X,Y,Z,double(BW(:,:,zStep(ii):zStep(ii+1))),Xq,Yq,Zq,'cubic');
    BW_interp =imbinarize(I_interp);
    if isempty(BW_stack)
        BW_stack =logical([BW_stack,BW_interp]);
    else
        BW_stack =cat(3,BW_stack,BW_interp);
    end
    size(BW_stack,3)
end
zC =round(length(zStep(end-1):zStep(end))*z/res); % Convert image to desired resolution
[X,Y,Z] = meshgrid(1:1:size(BW,2),1:1:size(BW,1),1:1:length(zStep(end-1):zStep(end)));
zL =linspace(1,length(zStep(end-1):zStep(end)),zC);
[Xq,Yq,Zq] = meshgrid(xL,yL,zL);
I_interp = interp3(X,Y,Z,double(BW(:,:,zStep(end-1):zStep(end))),Xq,Yq,Zq,'cubic');
BW_interp =imbinarize(I_interp);
BW_stack =cat(3,BW_stack,BW_interp);

end


function filtVoxels =filtEstimateFnc(radius,res,z)
% Estimate artefact size based on a defined radius of a perfect sphere.
% Sintax:
%     filtVoxels =filtEstimateFnc(radius,res,z)
% Inputs:
%     radius,           detected objects a radius length less than this value
%                       is removed
%     xyRes,            Resolusion in xy-axis of the images
%     z,                Resolution in the z-axis

% Outputs:
%     filtVoxels,       Defined noise size of objects with this number of
%                       voxels

% Filter noise less than 2.1904 ï¿½m^3 => 37 voxels default settings
% Default radius is 0.8 ï¿½m
if isempty(radius)
    radius =0.8; % decide a radius in ï¿½m
    SphereVol = 4/3*pi*radius^3;
    voxel = res/1000*res/1000*z/1000;
    filtVoxels = ceil(SphereVol./voxel);
else
    SphereVol = 4/3*pi*radius^3;
    voxel = res/1000*res/1000*z/1000;
    % estimate number of necessary voxels
    filtVoxels = round(SphereVol./voxel);
end

% SphereVol = 4/3*pi*r^3
% Estimate volume with a radius of 3ï¿½m
% SphereVol = 4/3*pi*3^3 = 113ï¿½m^3
% voxel of pixels = 0.272*0.272*0.8ï¿½m = 0.06ï¿½m^3/voxel
% Estimere antal voxel som er nï¿½dvendige med en radius af 3ï¿½m
% voxel = 113ï¿½m^3./0.06ï¿½m^3/voxel = 1884
% tr= (4/3*pi*7^3)/0.06;
% length(find(blobfilt.Volume>tr))
end

function [X,Y,Xq,Yq] =downScale2D(ImgGray,res,xy)
% Downsclae images
% Sintax:
%     [X,Y,Xq,Yq] =downScale2D(ImgGray,res,xy)
% Inputs:
%     ImgGray,      Grayscale image
%     res,          Output voxel resolution
%     xy,           xy pixel resolution of images

% Outputs:
%     X,            returns 2-D grid coordinates based on the coordinates
%                   contained in vectors x.
%     Y,            returns 2-D grid coordinates based on the coordinates
%                   contained in vectors y.
%     Xq,           returns 2-D grid coordinates based on the coordinates
%                   contained in vectors x.
%     Yq,           returns 2-D grid coordinates based on the coordinates
%                   contained in vectors y.

% dont change z
xC =round(size(ImgGray,2)*xy/res); % Convert 272 to 1000 nm pr pixel value
yC =round(size(ImgGray,1)*xy/res); % Convert 272 to 1000 nm pr pixel value

[X,Y] = meshgrid(1:1:size(ImgGray,2),1:1:size(ImgGray,1));
xL =linspace(1,size(ImgGray,2),xC);
yL =linspace(1,size(ImgGray,1),yC);

[Xq,Yq] = meshgrid(xL,yL);

end


function [I,Ic,cRect] = readImgFnc(folder1,folder2,folder3,fileName,ImgType,res,xy,buffer,CROP)
% Read stack of images
% Sintax:
%     [I,Ic,cRect] = readImgFnc(folder1,folder2,folder3,fileName,ImgType,res,xy,buffer,CROP)
% Inputs:
%     folder1,      First folder name
%     folder2,      Second folder name
%     folder3,      Third folder name
%     fileName,     Name of image that we read from
%     ImgType,      Image type (eg. TIF, png, JPEG)
%     res,          Output voxel resolution
%     xy,           xy pixel resolution of images
%     buffer,       starting point from analysis the images (eg. buffer 20
%                   => analysise from image 20)
%     CROP,         Crop area if necessary

% Outputs:
%     I,            Image stack
%     Ic,           Cropped Image stack
%     cRect,        Size and position of the crop rectangle in spatial coordinates

srcFiles = dir([fullfile(folder1,'/',folder2,'/',folder3),'/*.',ImgType]);
path=srcFiles.folder;
path=[path '\'];
N =size(srcFiles,1);
Img =imfill(imread([path,fileName,num2str(1+buffer),'.',ImgType]),'holes');
[~, ~, numberOfColorChannels] = size(Img);

% Change to desired image resolution
if res==xy
    if numberOfColorChannels > 1
        % Convert to grayscale
        ImgGray =rgb2gray(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),N);
        I(:,:,1) = imfill(imbinarize(ImgGray),'holes');
        %%
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N)])
            I_ori = (rgb2gray((imread([path,fileName,num2str(i),'.',ImgType]))));
            I(:,:,i) =imfill(imbinarize(I_ori),'holes');
        end
        clc
        disp([num2str(N),' binary images are loaded '])
    else
        ImgGray =(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),N);
        I(:,:,1) = imfill(ImgGray,'holes');
        
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N+buffer)])
            I_ori= (((imread([path,fileName,num2str(i),'.',ImgType]))));
            I(:,:,i-buffer) =imfill(I_ori,'holes');
        end
        clc
        disp([num2str(N),' binary images are loaded '])
    end
else
    if numberOfColorChannels > 1
        %%%%%%%%%%%%%%%%%%%%%%%%%% Downscale image
        ImgGray =rgb2gray(Img);
        [X,Y,Xq,Yq] =downScale2D(ImgGray,res,xy);
        Idown = interp2(X,Y,double(ImgGray),Xq,Yq);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(Idown,1),size(Idown,2),N);
        I(:,:,1) = imbinarize(Idown);
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N+buffer)])
            I_ori = (rgb2gray((imread([path,fileName,num2str(i),'.',ImgType]))));
            Idown = interp2(X,Y,double(I_ori),Xq,Yq);
            I(:,:,i-buffer) =imbinarize(Idown);
        end
        clc
        disp([num2str(N),' binary images are loaded '])
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%% Downscale image
        ImgGray =(Img);
        [X,Y,Xq,Yq] =downScale2D(ImgGray,res,xy);
        Idown = interp2(X,Y,double(ImgGray),Xq,Yq);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(Idown,1),size(Idown,2),N);
        I(:,:,1) = imbinarize(Idown);
        
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N+buffer)])
            I_ori= (((imread([path,fileName,num2str(i),'.',ImgType]))));
            Idown = interp2(X,Y,double(I_ori),Xq,Yq);
            I(:,:,i-buffer) =imbinarize(Idown);
        end
        clc
        disp([num2str(N),' binary images are loaded '])
    end
end
%%
% Crop Image
if CROP
    [Ic,cRect] = cropImg(I,CROP);
    
else
    Ic = [];
    cRect =[];
    
end

end

function BW = filtFillFnc(I, filtVoxels)
% Fill the holes of segmented images and remove small artefacts
% Sintax:
%     BW = filtFillFnc(I, filtVoxels)
% Inputs:
%     I,                Image stack
%     filtVoxels,       Defined noise size of objects with this number of
%                       voxels

% Outputs:
%     BW,               Image stack converted to binary form
BW = bwareaopen(I, filtVoxels);
% Fill the holes of segmented images
for i = 1:size(BW,3)
    BW(:,:,i) = imfill(BW(:,:,i),'holes');
end
end

function [x,cRect] = cropImg(I,CROP)
% Crop image stack
% Sintax:
%     [x,cRect] = cropImg(I,CROP)
% Inputs:
%     I,                Image stack
%     CROP,             Deside to crop or not

% Outputs:
%     x,                Cropped image stack
%     cRect,            Size and position of the crop rectangle in spatial coordinates

if CROP
    % Manuel Crop the image
    [~ ,cRect] = imcrop(imshow(I(:,:,23)));
    cRect = round(cRect);
    x0 =(imcrop(I(:,:,1),cRect));
    x = (zeros(size(x0,1),size(x0,2),size(I,3)));
    for i = 1:size(I,3)
        x(:,:,i) = imcrop(I(:,:,i),cRect);
    end
    close
    disp('Crop done')
else
    x = I;% Align the whole image and not just a subset
end

end
function [blobRemove,CCfilt] = blobFilt(blobinfo,CC,BW)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobRemove,CCfilt] = blobFilt(blobinfo,CC,BW)
% Inputs:
%     blobinfo,         Table with information of detected objects
%     CC,               Connected components, returned as a structure
%     BW,               Binary image stack

% Outputs:
%     blobRemove,       Filtered table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure
%% Remove top and bottom centroids
remv0 = blobinfo.Centroid(:,3)<3.5; % Exclude points from the first 3 image stacks
remv1 = blobinfo.Centroid(:,3)>(size(BW,3)-3.5); % Exclude points from the last 3 image stacks
remv2 =blobinfo.BoundingBox(:,6)<=1; % Remove dected objects in only a single plane
remv =logical(remv0+remv1+remv2); % Index for excluded points
CC.PixelIdxList(remv)=[] ;
blobinfo(remv,:) = [];
CC.NumObjects=size(blobinfo,1) ;
blobRemove = blobinfo;
CCfilt = CC;
end

function [blobfiltout,CCfilt] = blobParameters(blobRemove,CCfilt,z,res,side,window)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobfiltout,CCfilt] = blobParameters(blobRemove,CCfilt,z,res)
% Inputs:
%     blobRemove,       Filtered table with information of detected objects
%     CCfilt,           Filtered connected components, returned as a structure
%     z,                z resolution (distance between each image)
%     res,              Output voxel resolution
%     side,             Pial surface direction: left side 0, right side 1
%     window,           Defined window size of the imagestack

% Outputs:
%     blobfiltout,      Table with information of detected objects with
%                       added parameters
%     CCfilt,           Filtered Connected components, returned as a structure
%% Estimate the longest distance between centroid and all voxels
%%%%%% This value is not used in this study, but might be useful for other
%%%%%% studies
centroids =blobRemove.Centroid;
distMax =zeros(size(blobRemove,1),1);
diskLine =zeros(size(blobRemove,1),3);
for ii = 1:size(blobRemove,1)
    distArray = zeros(size(blobRemove.VoxelList{ii,1},1),1); %distance between centroid and voxels
    centroid1=centroids(ii,:);
    voxelList =blobRemove.VoxelList{ii,1};
    centroid1(3)=centroid1(3)*z/res;
    voxelList(:,3)=voxelList(:,3)*z/res;
    for j = 1:size(blobRemove.VoxelList{ii,1},1)
        distance_between_points = sqrt(sum((centroid1 - (voxelList(j,:))).^2)); % distance
        distArray(j) = distance_between_points;
    end
    % Find the farthest voxel from centroid
    distMax(ii)= max(distArray);
    idxMax = distArray==distMax(ii); % locical array of max distance
    distOne = find(idxMax, 1, 'first'); % if there more than one index with the length, then use the first one
    diskLine(ii,:)=voxelList(distOne,:);
end
%%%%% Estimate the distance between centroid and all voxels
blobRemove.DistMax = distMax*(res/1000);% Convert to micrometer
blobRemove.diskLine = diskLine;

%%%%% Estimate the average length of detected objects based on PrincipalAxisLength
adjustAxisLength =blobRemove.PrincipalAxisLength;
adjustAxisLength(:,3)=adjustAxisLength(:,3).*z/res; % Convert z to xy resolution pr pixel
%https://math.wikia.org/wiki/Ellipsoidal_quadratic_mean_radius
% Convert the length of axis to µm
adjustMicrometer =adjustAxisLength* (res/1000); % Convert length to micrometer
blobRemove.PrincipalAxisLengthAdjust = adjustMicrometer;
blobRemove.Diameter =sqrt(sum(blobRemove.PrincipalAxisLengthAdjust.^2,2)/3); % average/mean length of principalAxis
blobRemove.Radius = sqrt(sum((blobRemove.PrincipalAxisLengthAdjust/2).^2,2)/3); % average/mean length of principalAxis
%%%%% Estimate the width in the y and z plane 
width =sqrt(sum(blobRemove.PrincipalAxisLengthAdjust(:,2:3).^2,2)/2); % Width of objects (y and z-axis)
blobRemove.DistRatio =blobRemove.DistMax./blobRemove.Radius;
blobRemove.Width =width;

% Estimate Feret diameter
feretDiameter = zeros(size(blobRemove,1),1);
feretVolume = zeros(size(blobRemove,1),1);
feretPointMin = zeros(size(blobRemove,1),3);
feretPointMax = zeros(size(blobRemove,1),3);
for i = 1:size(blobRemove,1)
    try
        %distance between centroid and voxels
        voxelList =blobRemove.VoxelList{i,1};
        [k,fVol] = convhull(voxelList); % Get index of the convex position
        k = single(k);
        convexCoor =single([voxelList(k,1),voxelList(k,2),voxelList(k,3)]); % ConvexHull points
        convexCoor=unique(convexCoor,'rows'); % Remove redundant points
        convexCoor(:,3)=convexCoor(:,3)*z/1000;
        %         figure;
        %         scatter3(voxelList(:,1),voxelList(:,2),voxelList(:,3),'*')
        %         hold on
        %         plot3(voxelList(k,1),voxelList(k,2),voxelList(k,3))
        %         sqrt(sum((staticPoint-convexCoor).^2,2))
        %
        %         figure;
        %         scatter3(convexCoor(:,1),convexCoor(:,2),convexCoor(:,3))
        %         axis equal
        distArray = single(zeros(size(convexCoor,1)));
        for jj = 1:size(convexCoor,1)
            staticPoint =convexCoor(jj,:); %Define static point
            distArray(:,jj)=sqrt(sum((staticPoint-convexCoor).^2,2)); % Estimate distance of matrix
        end
        % Find the farthest voxel from centroid
        diskMax= max(max(distArray));
        idxMax = distArray==diskMax; % locical array of max distance
        distOne = find(idxMax, 1, 'first'); % if there more than one index with the length, then use the first one
        [r,c] = ind2sub(size(distArray),distOne);
        staticPointMin =convexCoor(c,:);
        lineMax =convexCoor(r,:);
        % % % % % Plot Line
        %         pts = [staticPointMin; lineMax];
        %         hold on
        %         % Because that's what line() wants to see
        %         plot3(pts(:,1), pts(:,2), pts(:,3),'green','LineWidth',3)
        %         axis equal
        
        if side==1
            %Depend on the side of the pial surface we should know the pointing
            %direction
            % Side ==1 pial surface to the right
            % Minimum side should be on the left side
            feretArray =[staticPointMin;lineMax];
            [~,idxferetMin]=min([feretArray(1);feretArray(2)]);
            [~,idxferetMax]=max([feretArray(1);feretArray(2)]);
            feretPointMin(i,:) =feretArray(idxferetMin,:); % save values
            feretPointMax(i,:)=feretArray(idxferetMax,:); % save values
        else
            % Side ==0 pial surface to the left
            % Minimum side should be on the right side where x-value is highest
            feretArray =[staticPointMin;lineMax];
            [~,idxferetMin]=min([feretArray(1);feretArray(2)]);
            [~,idxferetMax]=max([feretArray(1);feretArray(2)]);
            feretPointMin(i,:)=feretArray(idxferetMax,:); % save values
            feretPointMax(i,:) =feretArray(idxferetMin,:); % save values
        end
        feretDiameter(i,:)=diskMax; % save values
        feretVolume(i,:) =fVol;
        clear jj
    catch
        feretPointMin(i,:) =nan; % save values
        feretPointMax(i,:)=nan; % save values
        feretDiameter(i,:)=nan; % save values
        feretVolume(i,:) =nan;
        clear jj
    end
end

blobRemove.feretDiameter = feretDiameter*(res/1000);% Convert to micrometer
blobRemove.feretPointMin = feretPointMin;
blobRemove.feretPointMax = feretPointMax;
blobRemove.feretVolume = feretVolume;

%% Remove rows where feret diameter is empyty due to too few sections
remv = (isnan(feretPointMin(:,3))); %%%%%% here is the problem!!!!
blobRemove(remv,:) = [];
CCfilt.PixelIdxList(remv)=[] ;
CCfilt.NumObjects=size(blobRemove,1);
%% sort rows && Output
[blobfiltout, idxSort] = sortrows(blobRemove,1,'ascend');
%%%% Indext CCfilt
CCfilt.PixelIdxList=CCfilt.PixelIdxList(idxSort);
V=blobfiltout.Volume;
%%%% Estimate sphericity
blobfiltout.Sphericity = (pi^(1/3)*(6*V).^(2/3))./blobfiltout.SurfaceArea;%Sphericity
%%%% Estimate orientation
blobfiltout.thetaFeret =oriEstimateFeretFnc(window,blobfiltout,side);

end

function [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,blobfilt2D,CC2D,BW,kmeans3D,res,z,folder,name,SAVEFILES,SAVEIMAGES,SMOOTH)
% Visualize and save data of the 3D-reconstruction of pyramidal cells
% Sintax:
%     [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side,SAVEFILES,SMOOTH)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure
%     blobfilt2D        Filtered table with information of detected 2D-objects
%     CC2D              Filtered Connected components of 2D-images, returned as a structure
%     BW,               Binary image stack
%     kmeans3D,         kMean output data
%     res,              Output voxel resolution
%     z,                z pixel resolution of images
%     folder,           folder name
%     name,             name of subject/folder that will be analysed
%     SAVEFILES,        Save data (true or false)
%     SAVEIMAGES        Save images (true or false)
%     SMOOTH,           Smooth 3D-reconstruction images (true or false)

% Outputs:
%     pyramideNeuronImage,  Image stack of Pyramidal cells
%     roundNeuronImage,     Image stack of round neurons + artefacts

%% Define sampling window
window = [size(BW,2), size(BW,1),round(size(BW,3).*z/res)];  %[xLeft, yTop, length, width, hieght].
z_height = size(BW,3);
% clear BW;
roundNeurons=kmeans3D.RoundNeurons;
pyramidalNeurons=kmeans3D.PyramidalNeurons;
outlierNeurons=kmeans3D.outlierNeurons;
% Visualize neurons
neuronSmall = uint16(find(roundNeurons));
neuronBig =uint16(find(pyramidalNeurons));
blobfiltPyramid = blobfilt;
blobfiltPyramid(~(pyramidalNeurons),:) = [];

if SMOOTH
    labeledImage = labelmatrix(CCfilt);
    d = datestr(datetime('today'));
    % Only find the big and small volumes
    smallNeuronImage = ismember(labeledImage, neuronSmall);
    pyramideNeuronImage = ismember(labeledImage, neuronBig);
    
    % centroids
    cenAdRound = blobfilt.Centroid(neuronSmall,:); % finding the centriods for small neurons
    cenAdPyra = blobfilt.Centroid(neuronBig,:); % finding the centriods for pyramid neurons
    
    %%%%%%%%%%%% smooth the surface of the neurons for visualization
    dataRound = smooth3(smallNeuronImage,'box',3);
    dataPyra = smooth3(pyramideNeuronImage,'box',3);
    
    %%%%%% Visualize the 3D-reconstruction
    figure('units','normalized','outerposition',[0 0 1 1]);
    p = patch(isosurface((dataRound),0.5),...
        'FaceColor','interp','EdgeColor','interp');
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    alpha(0.45)                % set all patches transparency to 0.3
    daspect([1 1 1])
    view([-25,35])
    axis tight
    camlight
    lighting gouraud
    hold on;
    scatter3(cenAdRound(:,1),cenAdRound(:,2),cenAdRound(:,3),'filled','black')
    %%%%% Pyramidal cells
    hold on
    p = patch(isosurface((dataPyra),0.5),...
        'FaceColor','interp','EdgeColor','interp');
    axis tight
    p.FaceColor = 'green';
    p.EdgeColor = 'none';
    alpha(0.45)                % set all patches transparency to 0.3
    daspect([1 1 1])
    view([-25,35])
    axis tight
    camlight
    lighting gouraud
    
    hold on;
    % Plot Centroids
    scatter3(cenAdPyra(:,1),cenAdPyra(:,2),cenAdPyra(:,3),'filled','blue')
    
    %%%%%%%%%%% Plot Labels
    xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold'), zlabel('Z(µm)','FontSize',14,'FontWeight','bold')
    xticks(linspace(1,window(1),20));
    xticklabels(round(linspace(1,window(1),20)*res/1000));
    yticks(linspace(1,size(BW,1),20));
    yticklabels(round(linspace(1,window(2),20)*res/1000));
    zticks(linspace(1,size(BW,3),3));
    zticklabels(round(linspace(1,window(3),3)*res/1000));
    title(['Pyramid cells crop sequence ',d],'FontSize',14);
    set(gcf,'color','white');
       saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'m')
       saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'png')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 3D orientation plot
    % Visualize orientation of pyramidal cells
    for i =1:size(blobfiltPyramid,1)
        P1 = blobfiltPyramid.feretPointMin(i,:);
        P2 = blobfiltPyramid.feretPointMax(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','yellow', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    
        saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'m')
        saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'png')
else
    
    labeledImage = labelmatrix(CCfilt);
    d = datestr(datetime('today'));
    % Only find the big and small volumes
    smallNeuronImage = ismember(labeledImage, neuronSmall);
    pyramideNeuronImage = ismember(labeledImage, neuronBig);
    
    % centroids
    cenAdRound = blobfilt.Centroid(neuronSmall,:); % finding the centriods for small neurons
    cenAdPyra = blobfilt.Centroid(neuronBig,:); % finding the centriods for pyramid neurons
    
    %%%%%% Visualize the 3D-reconstruction
    figure;
    p = patch(isosurface((smallNeuronImage),0.5),...
        'FaceColor','interp','EdgeColor','interp');
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    alpha(0.45)                % set all patches transparency to 0.3
    daspect([1 1 1])
    view([-25,35])
    axis tight
    camlight
    lighting gouraud
    hold on;
    scatter3(cenAdRound(:,1),cenAdRound(:,2),cenAdRound(:,3),'filled','black')
    %%%%% Pyramidal cells
    hold on
    p = patch(isosurface((pyramideNeuronImage),0.5),...
        'FaceColor','interp','EdgeColor','interp');
    axis tight
    p.FaceColor = 'green';
    p.EdgeColor = 'none';
    alpha(0.45)                % set all patches transparency to 0.3
    daspect([1 1 1])
    view([-25,35])
    axis tight
    camlight
    lighting gouraud
    
    hold on;
    % Plot Centroids
    scatter3(cenAdPyra(:,1),cenAdPyra(:,2),cenAdPyra(:,3),'filled','blue')
    
    %%%%%%%%%%% Plot Labels
    xlabel('X(ï¿½m)','FontSize',14,'FontWeight','bold'), ylabel('Y(ï¿½m)','FontSize',14,'FontWeight','bold'), zlabel('Z(ï¿½m)','FontSize',14,'FontWeight','bold')
    xticks(linspace(1,size(labeledImage,2),8));
    xticklabels(round(linspace(1,window(1),8)*res/1000));
    yticks(linspace(1,size(labeledImage,1),8));
    yticklabels(round(linspace(1,window(2),8)*res/1000));
    zticks(linspace(1,size(labeledImage,3),3));
    zticklabels(round(linspace(1,window(3),3)*res/1000));
    title(['Pyramid cells crop sequence ',d],'FontSize',14);
    saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'png')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 3D orientation plot
    % Visualize orientation of pyramidal cells
    for i =1:size(blobfiltPyramid,1)
        P1 = blobfiltPyramid.feretPointMin(i,:);
        P2 = blobfiltPyramid.feretPointMax(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','yellow', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'png')
    close all
end
%%
if SAVEFILES==1
    %%
    %     Define the CCfilt map for pyramidal cells and non-pyramidal cells
    CCfiltAdjust =CCfilt;
    % remove objects that are not pyramidal cells
    CCfiltAdjust.PixelIdxList(~pyramidalNeurons)=[] ;
    CCfiltAdjust.NumObjects = size(CCfiltAdjust.PixelIdxList,2);
    CCfiltRound =CCfilt;
    % remove objects that are not pyramidal cells
    CCfiltRound.PixelIdxList(pyramidalNeurons)=[] ;
    CCfiltRound.NumObjects = size(CCfiltRound.PixelIdxList,2);
    %% Save files
    disp('Save variables')
    % kMeanOut = struct('Index',idx,'Clusters',C,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
    save(fullfile(folder, 'kmeans3D.mat'), 'kmeans3D')
    save(fullfile(folder, 'window.mat'), 'window')
    save(fullfile(folder, 'blobfilt.mat'), 'blobfilt', '-v7.3')
    save(fullfile(folder, 'blobfilt2D.mat'), 'blobfilt2D', '-v7.3')
    save(fullfile(folder, 'CCfilt.mat'), 'CCfilt')
    save(fullfile(folder, 'CCfiltAdjust.mat'), 'CCfiltAdjust')
    save(fullfile(folder, 'CCfiltRound.mat'), 'CCfiltRound')
    save(fullfile(folder, 'CC2D.mat'), 'CC2D')
    save(fullfile(folder, 'cenAdRound.mat'), 'cenAdRound')
    save(fullfile(folder, 'cenAdPyra.mat'), 'cenAdPyra')
    save(fullfile(folder, 'pyramidalNeurons.mat'), 'pyramidalNeurons')
    save(fullfile(folder, 'roundNeurons.mat'), 'roundNeurons')
    save(fullfile(folder, 'outlierNeurons.mat'), 'outlierNeurons')
    %%%%%%%%%% Save segmentation output
    save(fullfile(folder, 'blobfiltPyramid.mat'), 'blobfiltPyramid', '-v7.3')
    my_directory =fullfile(folder, ['centroidPyramidal_',name,'.xlsx']);
    xlswrite(my_directory,cenAdPyra); % Out of memory when you write
        clear kMeanOut blobfilt cenAdLow cenAdHigh blobfiltPyramid roundNeurons roundNeurons2
    pyramideNeuronImage =[];
    roundNeuronImage=[];
    %% Recontruct 3D image from labels bwconnected image
    if  SAVEIMAGES ==1
        disp('Reconstruct pyramide images')
        ImgSize =window(2)*window(1)*z_height;
        newImgPyra = false(ImgSize,1);
        for i = 1:size(CCfiltAdjust.PixelIdxList,2)
            sizeCC =size(CCfiltAdjust.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgPyra(CCfiltAdjust.PixelIdxList{1, i}(j)) =CCfiltAdjust.PixelIdxList{1, i}(j);
            end
        end
        pyramideNeuronImage = reshape(newImgPyra,[window(2),window(1),z_height]);
        save(fullfile(folder, 'pyramideNeuronImage.mat'), 'pyramideNeuronImage', '-v7.3')
        
        disp('Reconstruct Round images')
        ImgSize =window(2)*window(1)*z_height;
        newImgRound = false(ImgSize,1);
        for i = 1:size(CCfiltRound.PixelIdxList,2)
            sizeCC =size(CCfiltRound.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgRound(CCfiltRound.PixelIdxList{1, i}(j)) =CCfiltRound.PixelIdxList{1, i}(j);
            end
        end
        roundNeuronImage = reshape(newImgRound,[window(2),window(1),z_height]);
        save(fullfile(folder, 'roundNeuronImage.mat'), 'roundNeuronImage', '-v7.3')
    end
    
end
end


function [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,kmeansOutput,res,xy,z,folder,name,Dimension,SAVEIMAGES)
% Save data of the 3D-reconstruction of pyramidal cells
% Sintax:
%     [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,kmeansOutput,res,xy,z,folder,name,Dimension,SAVEIMAGES)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure
%     BW,               Binary image stack
%     kmeansOutput,     kMean output data
%     res,              Output voxel resolution
%     xy,               xy pixel resolution of images
%     z,                z pixel resolution of images
%     folder,           folder name
%     name,             name of subject/folder that will be analysed
%     Dimension,        Either 2D or 3D data  (2 or 3)
%     SAVEIMAGES        Save images (true/false)

% Outputs:
%     pyramideNeuronImage,  Image stack of Pyramidal cells
%     roundNeuronImage,     Image stack of round neurons + artefacts
%% Define sampling window
window = [size(BW,2), size(BW,1),round(size(BW,3).*z/res)];  %[xLeft, yTop, length, width, hieght].
z_height = size(BW,3);
clear BW;
%% Remove small neurons
pyramidalNeurons = kmeansOutput.PyramidalNeurons;
roundNeurons = kmeansOutput.RoundNeurons;
outlierNeurons=kmeansOutput.outlierNeurons;
blobfiltPyramid = blobfilt;
blobfiltPyramid(~(pyramidalNeurons),:) = [];
%% Estimate the two centroids for round and pyramidal cells
cenLow = blobfilt.Centroid(roundNeurons,:);
cenHigh = blobfilt.Centroid(pyramidalNeurons,:);
% Adjust centroids
cenAdLow = cenLow;
cenAdHigh = cenHigh;
cenAdHigh(:,3)=(cenAdHigh(:,3).*z/xy);% Convert Z to 272 nm pr pixel
cenAdLow(:,3)=cenAdLow(:,3).*z/xy;% Convert Z to 272 nm pr pixel
blobfiltPyramid.CentroidAdjust =cenAdHigh;
CCfiltAdjust =CCfilt;
% remove objects that are not pyramidal cells
CCfiltAdjust.PixelIdxList(~pyramidalNeurons)=[] ;
CCfiltAdjust.NumObjects = size(CCfiltAdjust.PixelIdxList,2);
CCfiltRound =CCfilt;
% remove objects that are not pyramidal cells
CCfiltRound.PixelIdxList(pyramidalNeurons)=[] ;
CCfiltRound.NumObjects = size(CCfiltRound.PixelIdxList,2);
%%
clear cenLow
clear cenHigh
if Dimension==3
    disp('Save variables')
    % kMeanOut = struct('Index',idx,'Clusters',C,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
    save(fullfile(folder, 'kmeans3D.mat'), 'kmeansOutput')
    save(fullfile(folder, 'window.mat'), 'window')
    save(fullfile(folder, 'CCfiltAdjust.mat'), 'CCfiltAdjust')
    save(fullfile(folder, 'CCfilt.mat'), 'CCfilt')
    save(fullfile(folder, 'CCfiltRound.mat'), 'CCfiltRound')
    save(fullfile(folder, 'blobfilt.mat'), 'blobfilt', '-v7.3')
    save(fullfile(folder, 'cenAdLow.mat'), 'cenAdLow')
    save(fullfile(folder, 'cenAdHigh.mat'), 'cenAdHigh')
    save(fullfile(folder, 'pyramidalNeurons.mat'), 'pyramidalNeurons')
    save(fullfile(folder, 'roundNeurons.mat'), 'roundNeurons')
    save(fullfile(folder, 'outlierNeurons.mat'), 'outlierNeurons')
    %%%%%%%%%% Save segmentation output
    save(fullfile(folder, 'blobfiltPyramid.mat'), 'blobfiltPyramid', '-v7.3')
    my_directory =fullfile(folder, ['centroidPyramidal_',name,'.xlsx']);
    xlswrite(my_directory,cenAdHigh); % Out of memory when you write
    clear CCfilt kMeanOut blobfilt cenAdLow cenAdHigh blobfiltPyramid roundNeurons roundNeurons2
    pyramideNeuronImage =[];
    roundNeuronImage=[];
    %% Recontruct 3D image from labels bwconnected image
    if  SAVEIMAGES ==1
        disp('Reconstruct pyramide images')
        ImgSize =window(2)*window(1)*z_height;
        newImgPyra = false(ImgSize,1);
        for i = 1:size(CCfiltAdjust.PixelIdxList,2)
            sizeCC =size(CCfiltAdjust.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgPyra(CCfiltAdjust.PixelIdxList{1, i}(j)) =CCfiltAdjust.PixelIdxList{1, i}(j);
            end
        end
        pyramideNeuronImage = reshape(newImgPyra,[window(2),window(1),z_height]);
        save(fullfile(folder, 'pyramideNeuronImage.mat'), 'pyramideNeuronImage', '-v7.3')
        
        disp('Reconstruct Round images')
        ImgSize =window(2)*window(1)*z_height;
        newImgRound = false(ImgSize,1);
        for i = 1:size(CCfiltRound.PixelIdxList,2)
            sizeCC =size(CCfiltRound.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgRound(CCfiltRound.PixelIdxList{1, i}(j)) =CCfiltRound.PixelIdxList{1, i}(j);
            end
        end
        roundNeuronImage = reshape(newImgRound,[window(2),window(1),z_height]);
        save(fullfile(folder, 'roundNeuronImage.mat'), 'roundNeuronImage', '-v7.3')
    end
else
    save(fullfile(folder, 'kmeans2D.mat'), 'kmeansOutput');
    blobfilt2D=blobfilt;
    save(fullfile(folder, 'blobfilt2D.mat'), 'blobfilt2D', '-v7.3')
    CCfilt2D =CCfilt;
    save(fullfile(folder, 'CCfilt2D.mat'), 'CCfilt2D')
    pyramideNeuronImage =[];
    roundNeuronImage=[];
    if  SAVEIMAGES ==1
        %% Recontruct 3D image from labels bwconnected image
        disp('Reconstruct pyramide images')
        ImgSize =window(2)*window(1)*z_height;
        newImgPyra = false(ImgSize,1);
        for i = 1:size(CCfiltAdjust.PixelIdxList,2)
            sizeCC =size(CCfiltAdjust.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgPyra(CCfiltAdjust.PixelIdxList{1, i}(j)) =CCfiltAdjust.PixelIdxList{1, i}(j);
            end
        end
        pyramideNeuronImage = reshape(newImgPyra,[window(2),window(1),z_height]);
        save(fullfile(folder, 'pyramideNeuronImage2D.mat'), 'pyramideNeuronImage', '-v7.3')
        
        disp('Reconstruct Round images')
        ImgSize =window(2)*window(1)*z_height;
        newImgRound = false(ImgSize,1);
        for i = 1:size(CCfiltRound.PixelIdxList,2)
            sizeCC =size(CCfiltRound.PixelIdxList{1, i},1);
            for j = 1:sizeCC
                newImgRound(CCfiltRound.PixelIdxList{1, i}(j)) =CCfiltRound.PixelIdxList{1, i}(j);
            end
        end
        roundNeuronImage = reshape(newImgRound,[window(2),window(1),z_height]);
        save(fullfile(folder, 'roundNeuronImage2D.mat'), 'roundNeuronImage', '-v7.3')
    end
end

end
function  saveData2DFnc(blobfiltProfile,kmeans2DP,folder)
% Save data of the 2D-data of pyramidal cells
% Sintax:
%     saveData2DFnc(blobfiltProfile,kmeans2DP,folder)
% Inputs:
%     blobfiltProfile,  Table with information of detected objects in 2D
%                       by using all cell profiles
%     kmeans2DP,         Struct of kmeans results
%     folder,           Folder location for all results
%% Remove small neurons
pyramidalNeurons = kmeans2DP.PyramidalNeurons;
roundNeurons = kmeans2DP.RoundNeurons;
blobfiltPyramid = blobfiltProfile;
blobfiltPyramid(~(pyramidalNeurons),:) = [];
blobfiltRound = blobfiltProfile;
blobfiltRound(~(roundNeurons),:) = [];
kmeans2DP.blobfiltProfile =blobfiltProfile;
kmeans2DP.blobfiltPyramid =blobfiltPyramid;
kmeans2DP.blobfiltRound =blobfiltRound;
%%
disp('Save variables')
% kMeanOut = struct('Index',idx,'Clusters',C,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
save(fullfile(folder, 'kmeans2DP.mat'), 'kmeans2DP')
save(fullfile(folder, 'blobfiltProfile.mat'), 'blobfiltProfile')
end


function theta =oriEstimateFeretFnc(window,blobfiltPyramid,side)
% Save data of the 3D-reconstruction of pyramidal cells
% Sintax:
%     [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,xy,z,folder,kMeanOut,name,side)
% Inputs:
%     window,           Size of the image stack window/frame
%     blobfiltPyramid,  Filtered table with information of detected
%     side,             Pial surface direction: left side 0, right side 1

% Outputs:
%     theta,            Measured orientation of cells in degrees

% side left = 0, right side 1
if side==0
    a = [window(1) 0 0];
    b = [1 0 0]; % Window(1) is the length of window
    u0 = (b-a)/(sqrt(sum(b-a).^2)); % u0 = [1 0 0]
    u0 =repmat(u0,size(blobfiltPyramid,1),1);
    
    
    c = blobfiltPyramid.feretPointMin; % position of centroid
    d =blobfiltPyramid.feretPointMax; % position of voxel farthest away from centroid
    dc = (d-c);
    dc_abs=sqrt(sum(dc.^2,2));
    u = dc./dc_abs;
    
    % Find angle between magnitude vector and all values
    dotU = dot(u,u0,2); % dot product
    uM =sqrt(sum(u.^2,2));
    u0M =sqrt(sum(u0.^2,2));
    theta= acos(dotU./(uM.*u0M))*180/pi;%angle 131.64
else
    a = [1 0 0];
    b = [window(1) 0 0]; % Window(1) is the length of window in x-direction
    u0 = (b-a)/(sqrt(sum(b-a).^2)); % u0 = [1 0 0]
    u0 =repmat(u0,size(blobfiltPyramid,1),1);
    c = blobfiltPyramid.feretPointMin; % position of centroid
    d =blobfiltPyramid.feretPointMax; % position of voxel farthest away from centroid
    dc = (d-c);
    dc_abs=sqrt(sum(dc.^2,2));
    u = dc./dc_abs;
    
    % Find angle between magnitude vector and all values
    dotU = dot(u,u0,2); % dot product
    uM =sqrt(sum(u.^2,2));
    u0M =sqrt(sum(u0.^2,2));
    theta= acos(dotU./(uM.*u0M))*180/pi;%angle 131.64
end
%         dotU = dot(u(1,:),u0); % -0.0414
%         uM =sqrt(sum(u(1,:).^2)); % 0.2200
%         u0M =sqrt(sum(u0.^2,2)); % 1
%         theta= acos(dotU./(uM.*u0M))*180/pi%angle 131.64
%% dot product
%         u = [3 -4 5]
%         v = [2 7 -3]
%
%         dotU =sum(u.*v)
%         uM =sqrt(sum(u.^2))
%         vM =sqrt(sum(v.^2))
%         theta= acos(dotU/(uM*vM))*180/pi%angle 131.64

end

function [CC2D,blobfilt2D] =Analysis2D2(BW,blobfilt,CCfilt,res,side)
%% 2D-analysis based on the 3D-reconstructed images
% Sintax:
%     [CC2D,blobfilt2D] =Analysis2D2(BW,blobfilt,CCfilt,res,side)
% Inputs:
%     BW,           Binary image stack
%     blobfilt,     Table with information of detected objects
%     CCfilt,       Filtered Connected components, returned as a structure
%     res,          Output voxel resolution
%     side,         Pial surface direction: left side 0, right side 1

% Outputs:
%     CC2D,         Filtered Connected components from 2D-sections
%     blobfilt2D,   Table with information of detected objects in 2D

% 2D analysis
CC2D = CCfilt;
disp('2D Analysis')
for i = 1:size(CC2D.PixelIdxList,2)
    % for i=1
    voxelList =blobfilt.VoxelList{i,1};
    % Detect number of detected slices that are in the same profile plane
    A = voxelList(:,3); % Only look at z-axis numbers
    groupC =groupcounts(A); %  Count number of equal number
    groupMax =max(groupC); % Find section with most numbers appears
    groupIdx =find(groupC==groupMax, 1, 'last' ); % Find the indices of where most sections appears
    [~, IdxArray] = unique(A, 'first'); % Index position of where they are placed
    if groupIdx==size(groupC,1)
        PixelIdxList2D =CC2D.PixelIdxList{1, i}(IdxArray(groupIdx):end,:);
    else
        PixelIdxList2D =CC2D.PixelIdxList{1, i}(IdxArray(groupIdx):IdxArray(groupIdx+1)-1,:);
    end
    % Compute the new 2D value
    CC2D.PixelIdxList{1, i} =PixelIdxList2D;
end
blobfilt2D =regionprops('table',CC2D,'Area','Centroid','BoundingBox','PixelList');
%% Feret estimation of a cell
feretDiameter = zeros(size(blobfilt2D,1),1);
feretPointMin = zeros(size(blobfilt2D,1),3);
feretPointMax = zeros(size(blobfilt2D,1),3);
SterelogyDistance = zeros(size(blobfilt2D,1),1);
SterelogyDetectLines = zeros(size(blobfilt2D,1),1);
distCenMin= zeros(size(blobfilt2D,1),1);
distCenMax= zeros(size(blobfilt2D,1),1);
distCenRatio= zeros(size(blobfilt2D,1),1);
c = clock;
disp(['2D processing have started. Time: ',num2str(c(4)),':',num2str(c(5))]);
%%
for i = 1:size(blobfilt2D,1)
    if i ==round(size(blobfilt2D,1)*0.1)
        c = clock;
        disp(['10% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.2)
        c = clock;
        disp(['20% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.3)
        c = clock;
        disp(['30% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.4)
        c = clock;
        disp(['40% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.5)
        c = clock;
        disp(['50% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.6)
        c = clock;
        disp(['60% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.7)
        c = clock;
        disp(['70% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.8)
        c = clock;
        disp(['80% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    elseif i ==round(size(blobfilt2D,1)*0.9)
        c = clock;
        disp(['90% Done. Time: ',num2str(c(4)),':',num2str(c(5))]);
    end
 
    try
        %Distance between centroid and voxels
        pixelList =blobfilt2D.PixelList{i,1};
        x=pixelList(:,1);
        y=pixelList(:,2);
        BinaryImg =false(size(BW,1),size(BW,2)); %% Generate image
        [ind] = sub2ind(size(BinaryImg),y,x); % get the indecies
        BinaryImg(ind)=1; % Place the indecies of the label image
        BinaryImg=imfill(BinaryImg,'holes');
        
        %%%%%%%%%% Define centroid %%% see if it is outside or inside the cell
        blob2D =regionprops('table',BinaryImg,'Area','Centroid','PixelList','Circularity','EquivDiameter');
        [~,Idx2D]=max(blob2D.EquivDiameter);
        % Remove other rows if there are any
        blob2D=blob2D(Idx2D,:);
        x=blob2D.PixelList{1,1}(:,1) ;
        y=blob2D.PixelList{1,1}(:,2);
        corners = [x y];
        
        BinaryImg =false(size(BW,1),size(BW,2)); %% Generate image
        [ind] = sub2ind(size(BinaryImg),y,x); % get the indecies
        BinaryImg(ind)=1; % Place the indecies of the label image
        %         imshow(BinaryImg)
        corners = permute(corners,[1 3 2]);
        corners = reshape(corners,[],2);
        corners = unique(corners,'rows');
        %%%%%%%%%%%%%%%%%%%%Plot
        %             plot(corners(:,1),corners(:,2),'sr','MarkerSize',5)
        %             hold on
        k = convhull(corners);
        hull_corners = corners(k,:);
        
        %%%%%%%%%%%%%%%%Plot
        %         plot(hull_corners(:,1),hull_corners(:,2),'r','LineWidth',3)
        %         plot(hull_corners(:,1),hull_corners(:,2),'ro','MarkerSize',10,'MarkerFaceColor','r')
        
        dx = hull_corners(:,1) - hull_corners(:,1)';
        dy = hull_corners(:,2) - hull_corners(:,2)';
        pairwise_dist = hypot(dx,dy);
        [max_dist,j] = max(pairwise_dist(:));
        [k1,k2] = ind2sub(size(pairwise_dist),j);
        
        point1 = hull_corners(k1,:);
        point2 = hull_corners(k2,:);
        %%%%%%%%%%%%%%%%Plot
        %         hold on
        %         plot([point1(1) point2(1)],[point1(2) point2(2)],'-db','LineWidth',3,'MarkerFaceColor','b')
        %         hold off
        staticPointMin = [point1,pixelList(1,3)];
        lineMax = [point2,pixelList(1,3)];
        if side==1
            %Depend on the side of the pial surface we should know the pointing
            %direction
            % Side ==1 pial surface to the right
            % Minimum side should be on the left side
            feretArray =[staticPointMin;lineMax];
            [~,idxferetMin]=min([feretArray(1);feretArray(2)]);
            [~,idxferetMax]=max([feretArray(1);feretArray(2)]);
            feretPointMin(i,:) =feretArray(idxferetMin,:); % save values
            feretPointMax(i,:)=feretArray(idxferetMax,:); % save values
            
            distCenMin(i) =pdist([blobfilt2D.Centroid(i,:);feretPointMin(i,:)],'euclidean');
            distCenMax(i) =pdist([blobfilt2D.Centroid(i,:);feretPointMax(i,:)],'euclidean');
            distCenRatio(i)= distCenMax(i)/distCenMin(i);
        else
            % Side ==0 pial surface to the left
            % Minimum side should be on the right side where x-value is highest
            feretArray =[staticPointMin;lineMax];
            [~,idxferetMin]=min([feretArray(1);feretArray(2)]);
            [~,idxferetMax]=max([feretArray(1);feretArray(2)]);
            % Opposite now
            % X-values closer to 0(min) are pointing toward the pial surface
            feretPointMin(i,:)=feretArray(idxferetMax,:); % save values
            feretPointMax(i,:) =feretArray(idxferetMin,:); % save values
            distCenMin(i) =pdist([blobfilt2D.Centroid(i,:);feretPointMin(i,:)],'euclidean');
            distCenMax(i) =pdist([blobfilt2D.Centroid(i,:);feretPointMax(i,:)],'euclidean');
            distCenRatio(i)= distCenMax(i)/distCenMin(i);
        end
        feretDiameter(i,:)=max_dist; % save values
        
        %%
        flag =1;
        count=0;
        LineNumber = 5; % Define the number of lines for the nucleator
        centroid1=blob2D.Centroid; % Centroid
        bigCellMean=0;
        %                     imshow(BinaryImg)
        
        while flag==1
            % %%%%%%%%% Generate image to get the boundary points
            dee =360/LineNumber; % Degree of rotation
            [yBoun, xBoun] = find(bwperim(BinaryImg)); %Find perimeter of objects in binary image
            % choose a random number from 1 to length
            idxRandom = randi(length(yBoun), 1);
            boundPoints = [xBoun,yBoun];
            distances = sqrt((centroid1(1)-boundPoints(:,1)).^2 + (centroid1(2)-boundPoints(:,2)).^2);
            [~,maxIdx]=(max(distances));
            
            disRatioX = xBoun(maxIdx)/xBoun(idxRandom);
            disRatioY = yBoun(maxIdx)/yBoun(idxRandom);
            %%%%%%%%%%%%%%%%%Plot
            %             plot(corners(:,1),corners(:,2),'sr','MarkerSize',5)
            %             axis equal
            %             hold on
            %             plot([centroid1(1) xBoun(idxRandom)],[centroid1(2) yBoun(idxRandom)],'-db','LineWidth',3,'MarkerFaceColor','b')
            %             scatter(xBoun(idxRandom),yBoun(idxRandom),'MarkerEdgeColor','black',...
            %                 'MarkerFaceColor','red',...
            %                 'LineWidth',1)
            %
            arrayPoints = zeros(LineNumber,2); % Array of boundary points for the nucleator
            arrayDistance = zeros(LineNumber,1); % Array of distance from points to centroid points
            arrayPoints(1,:)=[xBoun(idxRandom),yBoun(idxRandom)];
            arrayDistance(1,:) = (distances(idxRandom));        
            %%
            for ii = 2:LineNumber
                % Calculate rotation!!!!
                % aDegree = 30*pi/180;
                % ox = 0;
                % oy = 0;
                % px = 4;
                % py = 3;
                % qx = ox + cos(aDegree) * (px - ox) - sin(aDegree)*(py-oy)
                % qy = oy + sin(aDegree) * (px - ox) + cos(aDegree) * (py - oy)
                %%%%%%%%%%%%%% Always counter clockwise
                rotateAngle =dee*(ii-1);
                aDegree = (rotateAngle)*pi/180;
                xRatio = xBoun(idxRandom)*disRatioX;
                yRatio =yBoun(idxRandom)*disRatioY;
                %             qx = centroid1(1) + cos(aDegree) * (xBoun(idxRandom) - centroid1(1)) - sin(aDegree)*(yBoun(idxRandom)-centroid1(2));
                %             qy = centroid1(2) + sin(aDegree) * (xBoun(idxRandom) - centroid1(1)) + cos(aDegree) * (yBoun(idxRandom) - centroid1(2));
                qx = centroid1(1) + cos(aDegree) * (xRatio - centroid1(1)) - sin(aDegree)*(yRatio-centroid1(2));
                qy = centroid1(2) + sin(aDegree) * (xRatio - centroid1(1)) + cos(aDegree) * (yBoun(idxRandom) - centroid1(2));
                
                %%%%%%%%%%%%%%%%% plot
                %                 hold on
                %                 scatter(qx,qy,'b','filled')
                %                 plot([centroid1(1) qx],[centroid1(2) qy],'-db','LineWidth',3,'MarkerFaceColor','b')
                %                 axis equal

                %%%%%%%%%%%%% estimate slope
                a_slope =(qy-centroid1(2))/(qx-centroid1(1)); % estimate slope
                b =qy-a_slope*qx;
                data2=xBoun*a_slope+b; % Estimated y-data from the equation
                intersecPoint1 =find(floor(data2)==yBoun); % Find estimated y and compare to yBoun1
                intersecPoint2 =find(round(data2)==yBoun); % Find estimated y and compare to yBoun1
                intersecPoint3 =find(floor(data2)-1==yBoun); % Find estimated y and compare to yBoun1
                intersecPoint4 =find(floor(data2)+1==yBoun); % Find estimated y and compare to yBoun1
                intersecPoint =[intersecPoint1;intersecPoint2;intersecPoint3;intersecPoint4];
                intersecPoint=unique(intersecPoint);
                try
                    if isempty(intersecPoint)
                        intP0 =find(floor(data2)-3==yBoun); % Find estimated y and compare to yBoun1
                        intP1 =find(floor(data2)-2==yBoun); % Find estimated y and compare to yBoun1
                        intP2 =find(floor(data2)-1==yBoun); % Find estimated y and compare to yBoun1
                        intP3 =find(floor(data2)+1==yBoun); % Find estimated y and compare to yBoun1
                        intP4 =find(floor(data2)+2==yBoun); % Find estimated y and compare to yBoun1
                        intP5 =find(floor(data2)+3==yBoun); % Find estimated y and compare to yBoun1
                        intP6 =find(floor(data2)+4==yBoun); % Find estimated y and compare to yBoun1
                        intP7 =find(floor(data2)-4==yBoun); % Find estimated y and compare to yBoun1
                        intP8 =find(floor(data2)+5==yBoun); % Find estimated y and compare to yBoun1
                        intP9 =find(floor(data2)-5==yBoun); % Find estimated y and compare to yBoun1
                        intersecPoint =[intP1;intP2;intP3;intP4;intP5;intP6;intP7;intP8;intP9;intP0];
                        intersecPoint=unique(intersecPoint);
                        if qy>centroid1(2)%if the top point is above the centroid
                            idxCorrect =(yBoun(intersecPoint)>centroid1(2));
                            intersecPoint =intersecPoint(idxCorrect);
                            [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                            [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                            arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                            
                        else
                            %%
                            idxCorrect =(yBoun(intersecPoint)<centroid1(2));
                            intersecPoint =intersecPoint(idxCorrect);
                            [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                            [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                            arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                        end
                        
                    else
                        if qy>centroid1(2)%if the top point is above the centroid
                            %%
                            idxCorrect =(yBoun(intersecPoint)>centroid1(2));
                            intersecPoint =intersecPoint(idxCorrect);
                            [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                            %                 ThetaDiff =abs(rotateAngle-ThetaAngle); %estimate absolute difference
                            %                 [~,idxBound] =min(ThetaDiff);%the Find the angle with minimum difference
                            [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                            arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                        else
                            idxCorrect =(yBoun(intersecPoint)<centroid1(2));
                            intersecPoint =intersecPoint(idxCorrect);
                            [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                            %                 ThetaDiff =abs(rotateAngle-ThetaAngle); %estimate absolute difference
                            %                 [~,idxBound] =min(ThetaDiff);%the Find the angle with minimum difference
                            [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                            arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                            
                        end
                        
                    end
                    %%%%%%%%%%%%%%% plot intersection points position on
                    %%%%%%%%%%%%%%% boundary
                    %                     scatter(arrayPoints(ii,1),arrayPoints(ii,2),'MarkerEdgeColor','black',...
                    %                         'MarkerFaceColor','yellow',...
                    %                         'LineWidth',1)
                    %                     plot([centroid1(1) arrayPoints(ii,1)],[centroid1(2) arrayPoints(ii,2)],'b','LineWidth',3,'MarkerFaceColor','b')
                    % %
                    distances = sqrt((centroid1(1)-arrayPoints(ii,1)).^2 + (centroid1(2)-arrayPoints(ii,2)).^2);
                    arrayDistance(ii,:) = distances;
                catch
                    arrayDistance(ii,:) = 0;
                end
            end
            flag=0;  
            if count==0&&mean(arrayDistance)*(res/1000)>15&& length(nonzeros(arrayDistance))==LineNumber && bigCellMean==0
                bigCellMean = 1;
                arrayDistanceMatrix = zeros(5,5);
            end
            
            if bigCellMean==1
                flag =1;
                count = count+1;
                arrayDistanceMatrix(:,count)=arrayDistance;
                if count ==5
                    flag =0;
                    arrayDistanceMatrix(arrayDistanceMatrix == 0) = NaN;
                    arrayDistance =mean(arrayDistanceMatrix, 'omitnan')';% Set Zeros To â€˜NaNâ€™
                end
            end     
            if length(nonzeros(arrayDistance))<LineNumber && size(pixelList,1)>50 && bigCellMean==0
                flag =1;
                count = count+1;
                if count==5
                    LineNumber = 6;
                elseif count==10
                    if BinaryImg(round(centroid1(2)),round(centroid1(1)))==0% Check centroid position
                        % place centroid in the convex area
                        mPointRaw=((blob2D.PixelList{1,1}));
                        [N, edges]=histcounts(mPointRaw(:,1));
                        [~, idx]=max(N);
                        xCorr= ((edges(idx)+edges(idx+1))/2);
                        
                        [N, edges]=histcounts(mPointRaw(:,2));
                        [~, idx]=max(N);
                        yCorr= ((edges(idx)+edges(idx+1))/2);
                        
                        if BinaryImg(round(yCorr),round(xCorr))==1
                            centroid1(1)=xCorr;
                            centroid1(2)=yCorr;
                        end
                        LineNumber = 5;
                    end
                elseif count==15
                    LineNumber = 6;
                elseif count==26
                    flag=0;
                    if std(nonzeros(arrayDistance))>mean(nonzeros(arrayDistance))||size(find(arrayDistance==0),1)<5
                        BinaryPrim =bwperim(BinaryImg); % find the boundaries
                        blobFacts = regionprops(BinaryPrim, 'PixelIdxList'); % Find pixels at the boundary
                        row=[];
                        col =[];
                        if size(blobFacts,1)>1
                            
                            for jk = 1:size(blobFacts,1)
                                [Row,Col] = ind2sub([size(BinaryImg,1),size(BinaryImg,2)],blobFacts(jk).PixelIdxList);
                                row = [row;Row];
                                col = [col;Col];
                            end
                        else
                            [row,col] = ind2sub([size(BinaryImg,1),size(BinaryImg,2)],blobFacts.PixelIdxList);
                        end
                        cen=blob2D.Centroid;
                        distance_between_points = sqrt(sum((cen - ([col,row])).^2,2)); % distance
                        arrayDistance =mean(distance_between_points);     
                        if arrayDistance*(res/1000)>16.5 % if the radius is above 16.5µm then it is an outlier, as these objects are easy to detect
                            arrayDistance=nan;
                        end
                    end
                end
                
            end
            if count==10 &&mean(arrayDistance)*(res/1000)>15
                bigCellMean = 1;
                arrayDistanceMatrix = zeros(5,5);
                count = 0;
                flag =1;
            end
            
        end
        % Fill the table with values
        SterelogyDistance(i,:) = mean(nonzeros(arrayDistance));
        SterelogyDetectLines(i,:) = length(nonzeros(arrayDistance));
        blobfilt2D.Area(i) =blob2D.Area;
        blobfilt2D.Centroid(i,1:2) =blob2D.Centroid;
        blobfilt2D.Circularity(i) =blob2D.Circularity;
        blobfilt2D.PixelList(i) =blob2D.PixelList;
        blobfilt2D.AreaAdjust(i) =blobfilt2D.Area(i)*res/1000*res/1000;
        blobfilt2D.Equivdiameter(i) = sqrt(4*blobfilt2D.AreaAdjust(i)/pi);
        blobfilt2D.RadiusCircle(i) = sqrt(blobfilt2D.AreaAdjust(i)/pi); %A = pi*r2
    catch
        SterelogyDistance(i,:)=nan;
        SterelogyDetectLines(i,:)=nan;
        feretDiameter(i,:) =nan;
        feretPointMin(i,:)=nan;
        feretPointMax(i,:)=nan;
        distCenMin(i) =nan;
        distCenMax(i) =nan;
        distCenRatio(i)= nan; 
    end 
end
blobfilt2D.SterelogyDistanceRaw= SterelogyDistance;
blobfilt2D.SterelogyDistance= SterelogyDistance*(res/1000);
blobfilt2D.SterelogyDetectLines= SterelogyDetectLines;
blobfilt2D.feretDiameterRaw= feretDiameter;
blobfilt2D.feretDiameter= feretDiameter*(res/1000);
blobfilt2D.feretPointMin = feretPointMin;
blobfilt2D.feretPointMax = feretPointMax;
blobfilt2D.distCenMinRaw =distCenMin;
blobfilt2D.distCenMin =distCenMin*(res/1000);
blobfilt2D.distCenMaxRaw =distCenMax;
blobfilt2D.distCenMax =distCenMax*(res/1000);
blobfilt2D.distCenRatio= distCenRatio;
end

function [blobfiltProfile] =AnalysisAllProfiles2D(BW,blobfilt,CCfilt,res)
%% 2D analysis of All cell profiles
% Sintax:
%     [blobfiltProfile] =AnalysisAllProfiles2D(BW,blobfilt,CCfilt,res)
% Inputs:
%     BW,           Binary image stack
%     blobfilt,     Table with information of detected objects
%     CCfilt,       Filtered Connected components, returned as a structure
%     res,          Output voxel resolution

% Outputs:
%     blobfiltProfile,   Table with information of detected objects in 2D
%                        by using all cell profiles

blobfiltProfile = blobfilt;
CC2D = CCfilt;
blobfiltProfile(:,2:end)=[];
disp('2D Cell profile Analysis')
for i = 1:size(blobfilt,1)
    
    % for i=1:10
    % Finding connected components from CCfilt
    
    voxelList =blobfilt.VoxelList{i,1};
    % Detect number of detected slices that are in the same profile plane
    CC2D.NumObjects =1;
    CC2D.PixelIdxList =[];
    A = voxelList(:,3); % Only look at z-axis numbers
    groupC =groupcounts(A); %  Count number of equal number
    
    Nprofile =length(groupC);
    [~, IdxArray] = unique(A, 'first'); % Index position of where they are placed
    
    feretDiameter = zeros(Nprofile,1);
    profileSterelogyDistance = zeros(Nprofile,1);
    profileArea =zeros(Nprofile,1);
    profileEquivdiameter = zeros(Nprofile,1);
    for ij = 1:Nprofile
        %%
        try
            if ij==length(IdxArray)
                PixelIdxList2D =CCfilt.PixelIdxList{1, i}(IdxArray(ij):end,:);
            else
                PixelIdxList2D =CCfilt.PixelIdxList{1, i}(IdxArray(ij):IdxArray(ij+1),:);
            end
            CC2D.PixelIdxList{1} =PixelIdxList2D;
            bloblFilit2D =regionprops('table',CC2D,'Area','Centroid','BoundingBox','PixelList');
            
            %distance between centroid and voxels
            pixelList =bloblFilit2D.PixelList{1,1};
            
            x=pixelList(:,1);
            y=pixelList(:,2);
            BinaryImg =false(size(BW,1),size(BW,2)); %% Generate image
            [ind] = sub2ind(size(BinaryImg),y,x); % get the indecies
            BinaryImg(ind)=1; % Place the indecies of the label image
            %         imshow(BinaryImg)
            %         CC = bwconncomp(BW)
            
            blob2D =regionprops('table',BinaryImg,'Area','Centroid','PixelList','Circularity','EquivDiameter');
            [~,Idx2D]=max(blob2D.EquivDiameter);
            % Remove other rows if there are any
            blob2D=blob2D(Idx2D,:);
            x=blob2D.PixelList{1,1}(:,1) ;
            y=blob2D.PixelList{1,1}(:,2);
            corners = [x y];
            % Recreate image with only 1 profile
            BinaryImg =false(size(BW,1),size(BW,2)); %% Generate image
            [ind] = sub2ind(size(BinaryImg),y,x); % get the indecies
            BinaryImg(ind)=1; % Place the indecies of the label image
            %         imshow(BinaryImg)
            corners = permute(corners,[1 3 2]);
            corners = reshape(corners,[],2);
            corners = unique(corners,'rows');
            %%%%%%%%%%%%%%%%%%%%Plot
            %         plot(corners(:,1),corners(:,2),'sr','MarkerSize',5)
            %         hold on
            k = convhull(corners);
            hull_corners = corners(k,:);
            
            %%%%%%%%%%%%%%%%Plot
            %             plot(hull_corners(:,1),hull_corners(:,2),'r','LineWidth',3)
            %             plot(hull_corners(:,1),hull_corners(:,2),'ro','MarkerSize',10,'MarkerFaceColor','r')
            %
            dx = hull_corners(:,1) - hull_corners(:,1)';
            dy = hull_corners(:,2) - hull_corners(:,2)';
            pairwise_dist = hypot(dx,dy);
            [max_dist,~] = max(pairwise_dist(:));
            
            
            flag =1;
            count=0;
            while flag==1
                
                % %%%%%%%%% Generate image to get the boundary points
                LineNumber = 5;
                dee =360/LineNumber; % Degree of rotation
                
                %%%%%%%%%%%%%%%%% plot
                %                     plot(corners(:,1),corners(:,2),'sr','MarkerSize',5)
                %                         axis equal
                % hold on
                BinaryImg =false(size(BW,1),size(BW,2)); %% Generate image
                [ind] = sub2ind(size(BinaryImg),y,x); % get the indecies
                BinaryImg(ind)=1; % Place the indecies of the label image
                %             imshow(BinaryImg)
                [yBoun, xBoun] = find(bwperim(BinaryImg)); %Find perimeter of objects in binary image
                
                % choose a random number from 1 to length
                idxRandom = randi(length(yBoun), 1);
                % idxRandom = yBoun(idxRandom)
                boundPoints = [xBoun,yBoun];
                centroid1=blob2D.Centroid;
                distances = sqrt((centroid1(1)-boundPoints(:,1)).^2 + (centroid1(2)-boundPoints(:,2)).^2);
                % disRandom = distances(idxRandom);
                [~,maxIdx]=(max(distances));
                % idxRandom =maxIdx;
                
                disRatioX = xBoun(maxIdx)/xBoun(idxRandom);
                disRatioY = yBoun(maxIdx)/yBoun(idxRandom);
                %%%%%%%%%%%%%%%%%Plot
                %             hold on
                %             plot([centroid1(1) xBoun(idxRandom)],[centroid1(2) yBoun(idxRandom)],'-db','LineWidth',3,'MarkerFaceColor','b')
                %             scatter(xBoun(idxRandom),yBoun(idxRandom),'MarkerEdgeColor','black',...
                %                           'MarkerFaceColor','red',...
                %                           'LineWidth',1)
                
                arrayPoints = zeros(LineNumber,2);
                arrayDistance = zeros(LineNumber,1);
                arrayPoints(1,:)=[xBoun(idxRandom),yBoun(idxRandom)];
                arrayDistance(1,:) = (max(distances));
                
                %     %%%%%% Equation of the line chosen
                % a_slope =(yBoun1(idxRandom)-centroid1(2))/(xBoun1(idxRandom)-centroid1(1)); % estimate slope
                % b =yBoun1(idxRandom)-a_slope*xBoun1(idxRandom);
                % data2=xBoun1*a_slope+b; % Estimated y-data from the equation
                % intersecPoint =find(data2==yBoun1); % Find estimated y and compare to yBoun1
                %
                % interPoint = zeros(LineNumber,2);
                % interPoint(ii,:)=[xBoun1(intersecPoint),yBoun1(intersecPoint)]
                % intersectIdx=yBoun1==intersecPoint;
                %%
                for ii = 2:LineNumber
                    % Calculate rotation!!!!
                    % http://danceswithcode.net/engineeringnotes/rotations_in_2d/rotations_in_2d.html
                    % https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
                    % aDegree = 30*pi/180;
                    % ox = 0;
                    % oy = 0;
                    % px = 4;
                    % py = 3;
                    % qx = ox + cos(aDegree) * (px - ox) - sin(aDegree)*(py-oy)
                    % qy = oy + sin(aDegree) * (px - ox) + cos(aDegree) * (py - oy)
                    %%%%%%%%%%%%%% Always counter clockwise
                    rotateAngle =dee*(ii-1);
                    aDegree = (rotateAngle)*pi/180;
                    xRatio = xBoun(idxRandom)*disRatioX;
                    yRatio =yBoun(idxRandom)*disRatioY;
                    %  qx = centroid1(1) + cos(aDegree) * (xBoun(idxRandom) - centroid1(1)) - sin(aDegree)*(yBoun(idxRandom)-centroid1(2));
                    %  qy = centroid1(2) + sin(aDegree) * (xBoun(idxRandom) - centroid1(1)) + cos(aDegree) * (yBoun(idxRandom) - centroid1(2));
                    qx = centroid1(1) + cos(aDegree) * (xRatio - centroid1(1)) - sin(aDegree)*(yRatio-centroid1(2));
                    qy = centroid1(2) + sin(aDegree) * (xRatio - centroid1(1)) + cos(aDegree) * (yBoun(idxRandom) - centroid1(2));
                    
                    %%%%%%%%%%%%%%%%% plot
                    %             hold on
                    %             scatter(qx,qy,'b','filled')
                    %             plot([centroid1(1) qx],[centroid1(2) qy],'-db','LineWidth',3,'MarkerFaceColor','b')
                    %                         axis equal
                    
                    %%%%%%%%%%%%% estimate slope
                    a_slope =(qy-centroid1(2))/(qx-centroid1(1)); % estimate slope
                    b =qy-a_slope*qx;
                    data2=xBoun*a_slope+b; % Estimated y-data from the equation
                    intersecPoint1 =find(floor(data2)==yBoun); % Find estimated y and compare to yBoun1
                    intersecPoint2 =find(round(data2)==yBoun); % Find estimated y and compare to yBoun1
                    intersecPoint3 =find(floor(data2)-1==yBoun); % Find estimated y and compare to yBoun1
                    intersecPoint4 =find(floor(data2)+1==yBoun); % Find estimated y and compare to yBoun1
                    intersecPoint =[intersecPoint1;intersecPoint2;intersecPoint3;intersecPoint4];
                    %     intersecPoint =[intersecPoint1;intersecPoint2];
                    intersecPoint=unique(intersecPoint);
                    try
                        if isempty(intersecPoint)
                            intP0 =find(floor(data2)-3==yBoun); % Find estimated y and compare to yBoun1
                            intP1 =find(floor(data2)-2==yBoun); % Find estimated y and compare to yBoun1
                            intP2 =find(floor(data2)-1==yBoun); % Find estimated y and compare to yBoun1
                            intP3 =find(floor(data2)+1==yBoun); % Find estimated y and compare to yBoun1
                            intP4 =find(floor(data2)+2==yBoun); % Find estimated y and compare to yBoun1
                            intP5 =find(floor(data2)+3==yBoun); % Find estimated y and compare to yBoun1
                            intP6 =find(floor(data2)+4==yBoun); % Find estimated y and compare to yBoun1
                            intP7 =find(floor(data2)-4==yBoun); % Find estimated y and compare to yBoun1
                            intP8 =find(floor(data2)+5==yBoun); % Find estimated y and compare to yBoun1
                            intP9 =find(floor(data2)-5==yBoun); % Find estimated y and compare to yBoun1
                            intersecPoint =[intP1;intP2;intP3;intP4;intP5;intP6;intP7;intP8;intP9;intP0];
                            intersecPoint=unique(intersecPoint);
                            if qy>centroid1(2)%if the top point is above the centroid
                                idxCorrect =(yBoun(intersecPoint)>centroid1(2));
                                intersecPoint =intersecPoint(idxCorrect);
                                [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                                %                 ThetaDiff =abs(rotateAngle-ThetaAngle); %estimate absolute difference
                                %                 [~,idxBound] =min(ThetaDiff);%the Find the angle with minimum difference
                                [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                                arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                                
                            else
                                %%
                                idxCorrect =(yBoun(intersecPoint)<centroid1(2));
                                intersecPoint =intersecPoint(idxCorrect);
                                [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                                [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                                arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                                
                                
                                
                            end
                            
                        else
                            %%
                            if qy>centroid1(2)%if the top point is above the centroid
                                %%
                                idxCorrect =(yBoun(intersecPoint)>centroid1(2));
                                intersecPoint =intersecPoint(idxCorrect);
                                [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                                %                 ThetaDiff =abs(rotateAngle-ThetaAngle); %estimate absolute difference
                                %                 [~,idxBound] =min(ThetaDiff);%the Find the angle with minimum difference
                                [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                                arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                                
                                
                                
                            else
                                idxCorrect =(yBoun(intersecPoint)<centroid1(2));
                                intersecPoint =intersecPoint(idxCorrect);
                                [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints(1,:));
                                %                 ThetaDiff =abs(rotateAngle-ThetaAngle); %estimate absolute difference
                                %                 [~,idxBound] =min(ThetaDiff);%the Find the angle with minimum difference
                                [~,idxBound] =min(ThetaAngle/dee);%the Find the angle with minimum difference
                                arrayPoints(ii,:)=[xBoun(intersecPoint(idxBound)),yBoun(intersecPoint(idxBound))];
                                
                            end
                            
                        end
                        
                        %%%%%%%%%%%%%%%%% plot
                        %                             scatter(arrayPoints(ii,1),arrayPoints(ii,2),'MarkerEdgeColor','black',...
                        %                                 'MarkerFaceColor','yellow',...
                        %                                 'LineWidth',1)
                        %                             plot([centroid1(1) arrayPoints(ii,1)],[centroid1(2) arrayPoints(ii,2)],'b','LineWidth',3,'MarkerFaceColor','b')
                        %
                        distances = sqrt((centroid1(1)-arrayPoints(ii,1)).^2 + (centroid1(2)-arrayPoints(ii,2)).^2);
                        arrayDistance(ii,:) = distances;
                    catch
                        arrayDistance(ii,:) = 0;
                    end
                    
                    
                end
                flag=0;
                if length(nonzeros(arrayDistance))<LineNumber && size(pixelList,1)>50
                    flag =1;
                    count = [count+1];
                    if count>10
                        flag=0;
                    end
                end
                
            end
            feretDiameter(ij)=max_dist; % save values
            profileSterelogyDistance(ij,:) = mean(nonzeros(arrayDistance));
            profileArea(ij) =blob2D.Area;
            profileEquivdiameter(ij) = blob2D.EquivDiameter;
        catch
            feretDiameter(ij)=nan; % save values
            profileSterelogyDistance(ij,:) = nan;
            profileArea(ij) =nan;
            profileEquivdiameter(ij) = nan;
            
        end
        
    end
    profileData =struct('NumberProfiles',Nprofile,'profileFeretDiameter',feretDiameter,'profileSterelogyDistance',profileSterelogyDistance,...
        'profileArea',profileArea,'profileEquivdiameter',profileEquivdiameter);
    blobfiltProfile.Area(i) =mean(nonzeros(profileArea));
    blobfiltProfile.AreaAdjust(i) =blobfiltProfile.Area(i)*res/1000*res/1000;
    blobfiltProfile.Equivdiameter(i) = sqrt(4*blobfiltProfile.AreaAdjust(i)/pi);
    blobfiltProfile.RadiusCircle(i) = sqrt(blobfiltProfile.AreaAdjust(i)/pi); %A = pi*r2
    blobfiltProfile.SterelogyDistanceRaw(i)= mean(nonzeros(profileSterelogyDistance));
    blobfiltProfile.SterelogyDistance(i)= blobfiltProfile.SterelogyDistanceRaw(i)*(res/1000);
    blobfiltProfile.feretDiameterRaw(i)= mean(nonzeros(feretDiameter));
    blobfiltProfile.feretDiameter(i)= mean(nonzeros(feretDiameter))*(res/1000);
    blobfiltProfile.ProfileData(i) = profileData;
    
end
end







function [ThetaAngle] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints)
% Calculate the angle between two points
% Sintax:
%     [ThetaAngle,Phi] =calAngle(centroid1,intersecPoint,xBoun,yBoun,arrayPoints)
% Inputs:
%     centroid1,        Centroid coordinate
%     intersecPoint,    Point of intersection
%     xBoun,            x coordinate of boundary point
%     yBoun,            y coordinate of boundary point
%     arrayPoints,      Array of the boundary point of the nucleator

% Outputs:
%     ThetaAngle,       Measured angle in degrees

% Array of boundary points for the nucleator


% https://stackoverflow.com/questions/41518021/find-clockwise-angle-between-two-points-with-respect-to-an-arbitrary-origin
% A = ca and B = cb
% A = <ax - cx, ay - cy>
% B = <bx - cx, by - cy>
%
% Now that we have defined two vectors A & B lets find Theta.
% In order to find the cos(angle) between them we need to know both the magnitudes of A & B
% and the Dot Product between them.
% Magnitude or Length of a vector is sqrt( (i^2) + (j^2) )
% Dot Product is A*B = (Ai x Bi) + (Aj x Bj)
% Cos(Angle) =  (A*B) / (magA * magB)
% Theta = acos( Cos(Angle) )
% Phi = 360 - Theta

%%%- An example with points being a(5,7), b(3,5) & c(5,3)
% A = < 5-5, 7-3 > = <  0, 4 >
% B = < 3-5, 5-3 > = < -2, 2 >
%
% magA = sqrt( 0^2 + 4^2 ) = 4
% magB = sqrt( (-2)^2 + 2^2 ) = 2sqrt(2)
%
% cosAngle = A*B / (magA * magB)
% cosAngle = (0*-2 + 4*2) = 8 / ( 4 x 2sqrt(2) ) = 0.7071067812
%
% Theta = acos( cosAngle ) = 45 degrees
% Phi = 360 - 45 = 315

% A = arrayPoints(1,:);
% C =centroid1;
% a =[5,7];
% b =[3,5];
% c =[5,3];
%
% A =a-c;
% B =b-c;
% magA =sqrt(sum(A.^2));
% magB =sqrt(sum(B.^2));
% cosAngle = sum(A.*B)/ (magA * magB)
% ThetaAngle = acos(cosAngle)*180/pi;
% Phi = 360-ThetaAngle;

ThetaAngle=zeros(length(intersecPoint),1);
for k = 1:length(intersecPoint) 
    a =arrayPoints(1,:);
    b =[xBoun(intersecPoint(k)),yBoun(intersecPoint(k))];
    c =centroid1;
    
    A =a-c;
    B =b-c;
    magA =sqrt(sum(A.^2));
    magB =sqrt(sum(B.^2));
    cosAngle = sum(A.*B)/ (magA * magB);
    ThetaAngle(k) = acos(cosAngle)*180/pi;
end
end


function [kmeansOut] =kmeans2DFnc(bloblFilit2D,kmeans3D)
%% Values based on 2D-measurements
% Sintax:
%     [kmeansOut] =kmeans2DFnc(bloblFilit2D,kmeans3D)
% Inputs:
%     bloblFilit2D,  Filtered table with information of detected 2D-objects
%     kmeans3D,      Struct of kmeans results from 3D-measurements

% Outputs:
%     kmeansOut,     Struct of kmeans results

%% measure on the 2D-profile images
roundNeurons=kmeans3D.RoundNeurons;
pyramidalNeurons=kmeans3D.PyramidalNeurons;
outlierNeurons=kmeans3D.outlierNeurons;
idx =kmeans3D.idx;
C =kmeans3D.C;
kData =kmeans3D.kData;
Npyra = sum(pyramidalNeurons);
Nround = sum(roundNeurons);
Nout = sum(outlierNeurons);
Nratio = (size(bloblFilit2D,1)-Npyra)/size(bloblFilit2D,1);

kmeansOut = struct('idx',idx,'C',C,'kData',kData,'RoundNeurons',roundNeurons...
    ,'PyramidalNeurons',pyramidalNeurons,'outlierNeurons',outlierNeurons,...
    'NumberPyramid',Npyra,'NumberRound',Nround,'NumberOutliers',Nout,'NumberRatio',Nratio);
end

function [kmeansOut]=kmeans3DFnc(window,blobfilt,blobfilt2D,folder,res,z)
% Filter objects that are not pyramidal cells from the dataset
% Sintax:
%     [kmeansOut]=kmeans3DFnc(window,blobfilt,blobfilt2D,folder,res,z)
% Inputs:
%     window,           Window of image stack
%     blobfilt,     Table with information of detected objects
%     blobfilt2D,   Table with information of detected objects in 2D
%     folder,       Folder location for all results
%     res,          Output voxel resolution
%     z,            z resolution (distance between each image)

% Outputs:
%     kmeansOut,     Struct of kmeans results

% Define Outlier table
blobfilOutlier =blobfilt; % Find outliers
blobfilt2D.Diameter = blobfilt2D.SterelogyDistance*2; % Define diameter from 2D analysis
blobfilt2D.feretRatio =[blobfilt2D.feretDiameter]./[blobfilt2D.Diameter];
blobfilSmall2D =blobfilt2D;

% Remove outliers from Log normal distribution
posArtefacs = find(blobfilOutlier.feretDiameter); % Position of interesting points
posDiaArtefacs = find(blobfilOutlier.Diameter); % Position of interesting points
logNorFeret = log(blobfilOutlier.feretDiameter);
logFeret2D = log(blobfilSmall2D.feretDiameter);
logNorDia = log(blobfilSmall2D.Diameter);

%%%% Find outliers based on lognormal distribution
outCellsFeret = isoutlier(logNorFeret,'median'); %Detect outliers
outCellsFeret2D= isoutlier((logFeret2D),'median'); %Detect outliers
outCellsDia2D= isoutlier((logNorDia),'median'); %Detect outliers
% % 'median' (default) | 'mean' | 'quartiles' | 'grubbs' | 'gesd'

%%%%% Find only interesting points
removeOutlierFeret =posArtefacs(outCellsFeret); % Position of outliers in the table
removeOutlierFeret2D =posArtefacs(outCellsFeret2D); % Position of outliers in the table
removeOutlierDia2D =posDiaArtefacs(outCellsDia2D); % Position of outliers in the table

%%%% Remove redundant points
removeOutlier=unique([removeOutlierFeret;removeOutlierFeret2D;removeOutlierDia2D]);
outlierNeurons=false(size(blobfilt,1),1);
outlierNeurons(removeOutlier)=1;

%% %%%%%%%%% Find Non-pyramidal can pyramidal neurons
% Defining pyramidal cells
blobfilPyramid =blobfilt;
vol= nonzeros(blobfilPyramid.Volume*res/1000*res/1000*z/1000);
diaObj =(mean(vol));
posDia = find(vol<diaObj); % Position of interesting points
idxRadiusSmall=false(size(blobfilt,1),1); % Generete index
idxRadiusSmall(posDia)=1;
blobfilSmall= blobfilPyramid(idxRadiusSmall, :); % Table with only small radius values
% K means%%%%%%%%%%%%%%%%%%%%%%%%%%%
x= vol(idxRadiusSmall, :);
y= [blobfilSmall.Sphericity];
kData = [x y];
gm = fitgmdist(kData,2);
idx = cluster(gm,kData);
% cluster1 = (idx == 1); % |1| for cluster 1 membership
% cluster2 = (idx == 2); % |2| for cluster 2 membership
C =gm.mu;
roundC =C(:,1);
group =find(roundC==min(C(:,1)));
smallIAreaIdx =idx==group;
bigIdx=~smallIAreaIdx;

% Convert outliers into round group between values of detected centroid and max volume
thresMin =(min(C(:,1))+max(x(smallIAreaIdx)))/2;
idxRoundOutlier=x.*bigIdx<thresMin;
smallIAreaIdx=idxRoundOutlier;
bigIdx=~smallIAreaIdx;

idx = double(bigIdx*2);
idx =idx+smallIAreaIdx;
markerSize =4;
VISUALIZE =1;
if VISUALIZE
    figure, hold on
    lineWidth = 0.4;
    if group==1
        scatter(kData(idx==1,1),kData(idx==1,2),markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0.7 0 0.2],'MarkerFaceColor','r')
        scatter(kData(idx==2,1),kData(idx==2,2),markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0 0.5 1],'MarkerFaceColor','b')
    else
        scatter(kData(idx==1,1),kData(idx==1,2),markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0.7 0 0.2],'MarkerFaceColor','r')
        scatter(kData(idx==2,1),kData(idx==2,2),markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0 0.5 1],'MarkerFaceColor','b')
    end
    plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',2)
    xlabel('Volume (µm^3)'), ylabel('Sphericity')
    legend('Non-pyramidal','Pyramidal','Centroids',...
        'Location','NW')
    title 'Cluster Assignments and Centroids' 
end

set(gca,'FontSize',12);
saveas(gcf,fullfile(folder, 'kMeansplot3DVol'),'m')
saveas(gcf,fullfile(folder, 'kMeansplot3DVol'),'png')
set(0,'DefaultLegendAutoUpdate','off')
%%%%%% Add contour map
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
g = gca;
fcontour(gmPDF,[g.XLim g.YLim],'LineWidth',1.5)
plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',2)
saveas(gcf,fullfile(folder, 'kMeansplot3DVolContour'),'m')
saveas(gcf,fullfile(folder, 'kMeansplot3DVolContour'),'png')
hold off

%%%% Fine position of filtered cells and find their index
smallPos =posDia(smallIAreaIdx); % Position of the filteret neurons
nanValues=find(isnan(blobfilt2D.SterelogyDistance));
roundRemove = (unique([smallPos;nanValues])); % remove nan values as well
roundNeurons=false(size(blobfilt,1),1);
roundNeurons(roundRemove)=1;
roundNeurons = logical(roundNeurons);
pyramidalNeurons = ~(roundNeurons+outlierNeurons);

% Find intersection of values between round and outliwers
intersectArray = roundNeurons & outlierNeurons;
[rowsIntersect] = find(intersectArray);
% Remove values that are the same as filtered round neurons
outlierNeurons(rowsIntersect) =0;
dataRound= blobfilt(roundNeurons, :);
dataPyra= blobfilt(pyramidalNeurons, :);
dataOutlier= blobfilt(outlierNeurons, :);
dataRound2D= blobfilSmall2D(roundNeurons, :);
dataPyra2D= blobfilSmall2D(pyramidalNeurons, :);
dataOutlier2D= blobfilSmall2D(outlierNeurons, :);


Npyra = sum(pyramidalNeurons);
Nround = sum(roundNeurons);
Nout = sum(outlierNeurons);
Nratio = (size(blobfilSmall2D,1)-Npyra)/size(blobfilSmall2D,1);
% Plot kMeans figures include outliers
%# show points and clusters (color-coded)
markerSize =2;
figure, hold on
scatter(dataRound2D.Diameter,dataRound.feretDiameter,markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0.7 0 0.2],'MarkerFaceColor','r')
scatter(dataPyra2D.Diameter,dataPyra.feretDiameter,markerSize,'LineWidth',lineWidth,'MarkerEdgeColor',[0 0.5 1],'MarkerFaceColor','b')
scatter(dataOutlier2D.Diameter,dataOutlier.feretDiameter,15,'*','LineWidth',lineWidth,'MarkerEdgeColor','k')

xlabel('Diameter'), ylabel('FeretDiameter')
legend('Non-pyramidal','Pyramidal','Outliers',...
    'Location','NW')
title(['Cluster Assignments and Centroids. Ratio=',num2str(Nratio)])
hold off
saveas(gcf,fullfile(folder, 'kMeansplot'),'m')
saveas(gcf,fullfile(folder, 'kMeansplot'),'png')
hold off

%%%%%%%%%
den =round(sum(pyramidalNeurons)/(window(1)*window(2)*window(3))*10^9); % Estimate the density of neurons
kmeansOut = struct('idx',idx,'C',C,'kData',kData,'RoundNeurons',roundNeurons...
    ,'PyramidalNeurons',pyramidalNeurons,'outlierNeurons',outlierNeurons,...
    'NumberPyramid',Npyra,'NumberRound',Nround,'NumberOutliers',Nout,'NumberRatio',Nratio,'Window',window,'Density',den);
end


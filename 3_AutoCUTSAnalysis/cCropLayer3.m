% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read files %%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'Example';
folder1 ='AutuCUTS_Pipeline';
folder2 =[name,'_4_Prediction'];
folder3 =[name,'_UNetDense_stack'];
nameID = 'subject1';
fileName =[nameID,'pred'];
saveName = [nameID,'layer3'];
ImgType = 'png';
z = 800;     % z-hight between each image (nm)
xyRes = 272; % Resolution of images (nm)
numImage = [];
%% Run Program
% Create save image folders
[folder1,folder2,folder3]=fileInforFnc(folder1,folder2,folder3,ImgType);
[folder,folderData] =folderGenerateFnc(folder1,folder2,folder3);
%% Read images from the middle part of the tissue
[I,N] =readImgStackFnc(folder1,folder2,folder3,ImgType,fileName,numImage);
blobRemove =getCentroidsFnc(BW);
%% visOutput(blobRemove)
nBins = 15; % Define the binsize of the contour map
cRect = visDensityMapFnc(blobRemove,I,folderData,nBins,res);
%%%%%%%%%%%%%%%%%% Save images
for i = 1:N
    I_crop =readImgSingleLoopFnc(folder1,folder2,folder3,fileName,ImgType,cRect,i);
    imwrite(I_crop, [folder,'\',saveName,'_',num2str(i),'.png'],'png');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions
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
%     srcFiles,     Structure with file informations
%     folder1,      First folder name
%     folder2,      Second folder name
%     folder3,      Third folder name
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

function [folder,folderData] =folderGenerateFnc(folder1,folder2,folder3)
% Create the save folder of images
% Sintax:
%     [folder,folderData] =folderGenerateFnc(folder1,folder2,folder3)
% Inputs:
%     folder1,          First folder name
%     folder2,          Second folder name
%     folder3,          Third folder name

% Outputs:
%     folder,           Directory for folder of cropped images
%     folderData,       Directory for folder of cropped images data

% Create image folder
s = what(folder1);
savePath=[s.path,'\',folder2];
folderResutls =[folder3,'_layer3'];
folder = [savePath,'\',folderResutls];
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end
% Create data of layer 3 folder
s = what(folder1);
savePathData=[s.path,'\',folder2];
folderResutls =[folder3,'_layer3'];
folderData = [savePathData,'\',folderResutls,'_data'];
if ~exist(folderData, 'dir')
    mkdir(folderData);
end

end

function [I,N] =readImgStackFnc(folder1,folder2,folder3,ImgType,fileName,numImage)
% Read stack of images used to crop images
% Sintax:
%     [I,N] =readImg(folder1,folder2,folder3,ImgType,fileName,numImage)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     folder3,     Second folder name
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     fileName,    Name of image that we read from
%     numImage,    Set the number of images from the folder to be loaded

% Outputs:
%     I,            Output image stack
%     N,            Number of images in the folder

% Getting directory of images
srcFiles = dir([fullfile(folder1,'/',folder2,'/',folder3),'/*.',ImgType]);
%%
path=srcFiles.folder;
path=[path '\'];
N =size(srcFiles,1); % Number of total images

% Determine how many images we are going to read
% If we have not defined a number of images that should be read => read whole stack in folder
if isempty(numImage)
    numImage=N;
    flag = 1;
else
    numImage=round(N/3);
    flag = 0;
end

Img =imread([path,fileName,num2str(1),'.',ImgType]); % Read image number 1
[~, ~, numberOfColorChannels] = size(Img); % Detect if images are RBG or grayscale
if numberOfColorChannels > 1
%%%%%%%%%%%%%% RGB images
    if flag==1
        ImgGray =rgb2gray(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),numImage);
        I(:,:,1) = imbinarize(ImgGray);
        for i = 2:numImage
            disp(['Load image number ',num2str((i)),' of ',num2str((numImage)),' images',])
            I_ori = uint8(rgb2gray((imread([path,fileName,num2str((i)),'.',ImgType]))));
            I(:,:,i) =imbinarize(I_ori);
        end
        clc
        disp([numImage,' binary images are loaded '])
        
    else
        ImgGray =rgb2gray(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),numImage);
        I(:,:,1) = imbinarize(ImgGray);
        list = [1:round(numImage*0.33),round((N/2)*0.90):round((N/2)*0.90)+round(numImage*0.33), N-round(numImage*0.33):N];
        
        for i = 2:length(list)
            disp(['Load image number ',num2str(list(i)),'(',num2str((i)),') of ',num2str(length(list)),' images',])
            I_ori = uint8(rgb2gray((imread([path,fileName,num2str(list(i)),'.',ImgType]))));
            I(:,:,i) =imbinarize(I_ori);
        end
        clc
        disp([length(list),' binary images are loaded '])
        
    end
    
else
%%%%%%%%%%%%%% Grayscale images
    if flag==1
        ImgGray =rgb2gray(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),numImage);
        I(:,:,1) = imbinarize(ImgGray);
        for i = 2:numImage
            disp(['Load image number ',num2str((i)),' of ',num2str((numImage)),' images',])
            I_ori = uint8(rgb2gray((imread([path,fileName,num2str((i)),'.',ImgType]))));
            I(:,:,i) =imbinarize(I_ori);
        end
        clc
        disp([numImage,' binary images are loaded '])
    else
        ImgGray =logical(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),numImage);
        I(:,:,1) = ImgGray;
        list = [1:round(numImage*0.33),round((N/2)*0.90):round((N/2)*0.90)+round(numImage*0.33), N-round(numImage*0.33):N];
        
        for i = 2:length(list)
            disp(['Load image number ',num2str(list(i)),' of ',num2str(length(list)),' images',])
            I_ori = uint8(((imread([path,fileName,num2str(list(i)),'.',ImgType]))));
            I(:,:,i) =logical(I_ori);
        end
    end
end
clc
disp([num2str(numImage),' binary images are loaded '])
end



function blobRemove =getCentroidsFnc(BW)
%%
% Remove noise from images and get centroid points
% Sintax:
%     blobRemove =getCentroidsFnc(BW)
% Inputs:
%     BW,               Binary image stack
% Outputs:
%     blobRemove,       Table with information of detected objects

disp('Detect connections')
CC = bwconncomp(BW);
blobinfo = regionprops3(CC,'BoundingBox','Centroid','VoxelList');
[blobRemove, ~]= blobFilt(blobinfo,CC,BW);
disp('Connections loaded')
end
function cRect = visDensityMapFnc(blobRemove,I,folderData,nBins,res)
% Visualize the density plot of centroid placed in the stack of images
% Sintax:
%     cRect = visDensityMapFnc(blobRemove,I,folderData,cropLines,nBins)
% Inputs:
%     blobRemove,       blobRemove
%     I,                Image stack
%     folderData,       Directory for folder of cropped images data
%     nBins,            Desize the bin size of the density map
%     res,              Pixel resolution

% Outputs:
%     cRect,       Size and position of the crop rectangle in spatial coordinates
BW = I;
centroids =blobRemove.Centroid;

centroids(:,3)=centroids(:,3);
numberLabel = 8;
figure(2);
scatter3(centroids(:,1),centroids(:,2),centroids(:,3),'filled','blue')
title(['Centre of mass of ',num2str(size(centroids,1)),' filtered pyramidal cells.'],'FontSize',14);
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold'), zlabel('Z(µm)','FontSize',14,'FontWeight','bold')
xticks(round(linspace(1,size(BW,2),numberLabel)));
xticklabels(round(round(linspace(1,size(BW,2),numberLabel))*res/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*res/1000));
%% Assign axis values from centroids
X=centroids(:,1);
Y=centroids(:,2);
%% Interpolated image
minData = min([X;Y]);
maxData = max([X;Y]);
bins = linspace(minData, maxData, nBins);
N = histcounts2(Y(:), X(:), bins, bins);

subplot(1,2,1)
scatter(X, Y, 'b.');
axis equal;
xlim([1,size(BW,2)])
ylim([1,size(BW,1)])
title('xy-plane','FontSize',17,'FontWeight','bold')
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold')
xticks(round(linspace(1,size(BW,2),numberLabel)));
xticklabels(round(round(linspace(1,size(BW,2),numberLabel))*res/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*res/1000));

% Plot heatmap:
sh = subplot(1,2,2);
r = linspace(1,maxData,size(N,1));
pcolor(r,r,N);
shading interp
axis equal;
xlim([1,size(BW,2)])
ylim([1,size(BW,1)])
set(sh, 'YDir', 'normal')
title('xy-plane','FontSize',17,'FontWeight','bold')
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold')
xticks(linspace(1,size(BW,2),numberLabel));
xticklabels(round(linspace(1,size(BW,2),numberLabel)*res/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*res/1000));
%%% Save files
saveas(gcf,fullfile(folderData, 'DensityPlot'),'m')
saveas(gcf,fullfile(folderData, 'DensityPlot'),'pdf')
saveas(gcf,fullfile(folderData, 'DensityPlot'),'png')
%%%% Mark area

%% Visualize all datapoints
figure(3)
pcolor(r,r,N);
shading interp
xlim([1,size(BW,2)])
ylim([1,size(BW,1)])
xticks(linspace(1,size(BW,2),numberLabel));
xticklabels(round(linspace(1,size(BW,2),numberLabel)*res/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*res/1000));
title('xy-plane','FontSize',17,'FontWeight','bold')
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold')

title('Choose top left cornor','FontSize',17,'FontWeight','bold')
[xi, yi] = ginput(1);
% Visualize crop output
% p1 = [xi; yi];
% p2 = [xi; yi+height];
% p3 = [xi+width; yi];
% p4 = [xi+width; yi+height];
hold on
xi = round(xi);
y = yi;
lineWidth =2;
scatter(xi,y,'r','filled')
plot([xi xi],[0 yi],'--r','LineWidth',lineWidth)
title('Choose bottom right cornor','FontSize',17,'FontWeight','bold')
[x2, y2] = ginput(1);
% Visualize crop output
% p1 = [xi; yi];
% p2 = [xi; yi+height];
% p3 = [xi+width; yi];
% p4 = [xi+width; yi+height];
hold on
x2 = round(x2);
scatter(xi,y2,'r','filled') % bottom left
scatter(x2,y2,'r','filled') % bottom right
scatter(x2,y,'r','filled') % top right
plot([x2 x2],[y2 y],'--r','LineWidth',lineWidth) % right line
plot([xi x2],[yi yi],'--r','LineWidth',lineWidth) % bottom line
plot([xi x2],[y2 y2],'--r','LineWidth',lineWidth) % top line
cRect=round([xi y2 x2-xi y-y2]); % ROI of the crop image  [xmin ymin width height].
if sum(mod(cRect,2)==1)>0
    idx = find(mod(cRect,2)==1);
    cRect(idx)=cRect(idx)+1;
end
title('Region of Interest is within this square','FontSize',17,'FontWeight','bold')
%%
saveas(gcf,fullfile(folderData, 'cropPlot'),'m')
saveas(gcf,fullfile(folderData, 'cropPlot'),'pdf')
saveas(gcf,fullfile(folderData, 'cropPlot'),'png')
save(fullfile(folderData, 'blobRemove.mat'), 'blobRemove')
save(fullfile(folderData, 'cRect.mat'), 'cRect')

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
remv0 = blobinfo.Centroid(:,3)<3.5;
remv1 = blobinfo.Centroid(:,3)>(size(BW,3)-3.5);
% remv2 =blobinfo.BoundingBox(:,6)<=1;
remv =logical(remv0+remv1);
CC.PixelIdxList(remv)=[] ;
blobinfo(remv,:) = [];
CC.NumObjects=size(blobinfo,1) ;
blobRemove = blobinfo;
CCfilt = CC;
end

function [x] = readImgSingleLoopFnc(folder1,folder2,folder3,fileName,ImgType,cRect,i)
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
%     cRect,       Size and position of the crop rectangle in spatial coordinates
%     i,           Image number from loop

% Outputs:
%     x,                    Output image name

% Get directory
srcFiles = dir([fullfile(folder1,'/',folder2,'/',folder3),'/*.',ImgType]);
path=srcFiles.folder;
% Detect number of images
N2 = size(srcFiles,1);
path=[path '\'];
disp(['Load image number ',num2str(i),' out of ',num2str(N2)])
try
    I_ori = uint8(rgb2gray((imread([path,fileName,num2str(i),'.',ImgType]))));
    I=imbinarize(I_ori);
catch
    I_ori = uint8(((imread([path,fileName,num2str(i),'.',ImgType]))));
    I=imbinarize(I_ori);
end
% I = (imread([path,srcFiles(i).name]));

%% Crop Image
x = cropImgLoopFnc(I,cRect);
x = imfill(x,'holes');
end

function x = cropImgLoopFnc(I,cRect)
% Manuel Crop the image
cRect = round(cRect);
%     x = zeros(cRect(4)+1,cRect(3)+1,size(I,3),'uint8');

x = imcrop(I,cRect);

%% Remove rows and columns
if mod(size(x,1),2) ==1
    x(end,:,:)=[];
end
if mod(size(x,2),2) ==1
    x(:,end,:) =[];
end
end

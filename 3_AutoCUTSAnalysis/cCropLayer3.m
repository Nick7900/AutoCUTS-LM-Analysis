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
z = 900;     % z-hight between each image (nm)
xyRes = 272; % Resolution of images (nm)

%% Run Program
% Create save image folder
[folder,folderData] =folderGenerateFnc(folder1,folder2,folder3);
% Read images from the middle part of the tissue
[I,N] =readImgStackFnc(folder1,folder2,folder3,ImgType,fileName);
%%
% Remove noise and plot centroids
% filter small cells away
radius = 0.5; % estimate number of pixels required to reach a radius of the determine size.
filtVoxels =filtEstimateFnc(radius,xyRes,z); % remove objects with less pixels than filtVoxels
blobRemove =getCentroidsFnc(I,filtVoxels); % Measure properties of 3D volumetric objects

%% visOutput(blobRemove)
cropLines = 8; % Define number of lines that is used to crop the image stack
nBins = 15;
cRect = visDensityMapFnc(blobRemove,I,folderData,cropLines,nBins); % visualize density map

% Save images
for i = 1:N
    I_crop =readImgSingleLoopFnc(folder1,folder2,folder3,fileName,ImgType,cRect,i);
    imwrite(I_crop, [folder,'\',saveName,'_',num2str(i),'.png'],'png');
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,N] =readImgStackFnc(folder1,folder2,folder3,ImgType,fileName)
% Read stack of images used to crop images
% Sintax:
%     [I,N] =readImg(folder1,folder2,folder3,ImgType,fileName)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     folder3,     Second folder name
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     fileName,    Name of image that we read from

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
% If we read less than 200 images than read all images in folder
if N<200==1
    numImage=N;
    flag = 1;
else
    numImage=round(N/3);
    flag = 0;
end

Img =imread([path,fileName,num2str(1),'.',ImgType]); % Read image number 1
[~, ~, numberOfColorChannels] = size(Img);
%%
if numberOfColorChannels > 1
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
    %%%%%%%%%%%%%%%%%%%%%%%%%% Downscale image
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

function blobRemove =getCentroidsFnc(I,filtVoxels)
% Remove noise from images and get centroid points
% Sintax:
%     [I,N] =readImg(folder1,folder2,folder3,ImgType,fileName)
% Inputs:
%     I,                Image stack
%     filtVoxels,       Filter noise away less than this value

% Outputs:
%     blobRemove,       Table with information of detected objects

BW = bwareaopen(I, filtVoxels);
% Fill the holes of segmented images
for i = 1:size(BW,3)
    BW(:,:,i) = imfill(BW(:,:,i),'holes');
end
% Detect connected neurons
disp('Detect connections')
CC = bwconncomp(BW);
% Measure properties of 3D volumetric objects
blobinfo = regionprops3(CC,'BoundingBox','Centroid','Volume','Orientation','VoxelList','PrincipalAxisLength');

% Filter small objects away
[blobRemove, ~]= blobFiltFnc(blobinfo,CC,BW);
disp('Connections loaded')
end

function cRect = visDensityMapFnc(blobRemove,I,folderData,cropLines,nBins)
% Visualize the density plot of centroid placed in the stack of images
% Sintax:
%     cRect = visDensityMapFnc(blobRemove,I,folderData,cropLines,nBins)
% Inputs:
%     I,                Image stack
%     folderData,       Directory for folder of cropped images data
%     cropLines,        Define number of lines that need to appear
%     nBins,            Desize the bin size of the density map

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
xticklabels(round(round(linspace(1,size(BW,2),numberLabel))*272/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*272/1000));
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
xticklabels(round(round(linspace(1,size(BW,2),numberLabel))*272/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*272/1000));

% Plot heatmap:
sh = subplot(1,2,2);
r = linspace(1,maxData,size(N,1));
pcolor(r,r,N);
shading interp
%  colorbar
axis equal;
xlim([1,size(BW,2)])
ylim([1,size(BW,1)])
set(sh, 'YDir', 'normal')
title('xy-plane','FontSize',17,'FontWeight','bold')
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold')
xticks(linspace(1,size(BW,2),numberLabel));
xticklabels(round(linspace(1,size(BW,2),numberLabel)*272/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*272/1000));
%%% Save files
saveas(gcf,fullfile(folderData, 'DensityPlot'),'m')
saveas(gcf,fullfile(folderData, 'DensityPlot'),'png')
%% Visualize all datapoints
figure(3)
pcolor(r,r,N);
shading interp
xlim([1,size(BW,2)])
ylim([1,size(BW,1)])
xticks(linspace(1,size(BW,2),numberLabel));
xticklabels(round(linspace(1,size(BW,2),numberLabel)*272/1000));
yticks(linspace(1,size(BW,1),numberLabel));
yticklabels(round(linspace(1,size(BW,1),numberLabel)*272/1000));
title('xy-plane','FontSize',17,'FontWeight','bold')
xlabel('X(µm)','FontSize',14,'FontWeight','bold'), ylabel('Y(µm)','FontSize',14,'FontWeight','bold')

[xi, yi] = ginput(1);
% Visualize crop output
% p1 = [xi; yi];
% p2 = [xi; yi+height];
% p3 = [xi+width; yi];
% p4 = [xi+width; yi+height];
hold on
xi = round(xi);
y = yi;
scatter(xi,y,'r','filled')
plot([xi xi],[0 yi],'--r')

x2 = zeros(1,cropLines);
y2 = zeros(1,cropLines);
for i = 1:cropLines
    x2(i) = xi+256*(i*1+4) ;
    y2(i) = yi;
    plot([x2(i) x2(i)],[0 yi],'--b')
    str = num2str(i);
    text(x2(i)+100,y2(i),str)
end
scatter(x2,y2,'b','filled')
% x2 = [xi+256*4 xi+256*6 xi+256*8 xi+256*10];
% y2 = [yi yi yi yi];

%  point number of marked width
prompt = 'Where to crop the image from?\n ';
numMark = input(prompt);
cRect=[xi 0 x2(numMark)-xi size(BW,2)]; % ROI of the crop image
title(['Crop from line ',num2str(numMark)],'FontSize',17,'FontWeight','bold')

saveas(gcf,fullfile(folderData, 'cropPlot'),'m')
saveas(gcf,fullfile(folderData, 'cropPlot'),'png')
save(fullfile(folderData, 'blobRemove.mat'), 'blobRemove')
save(fullfile(folderData, 'cRect.mat'), 'cRect')
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

function [blobRemove,CCfilt] = blobFiltFnc(blobinfo,CC,BW)
% Filter the data of detected objects from the table.
% Sintax:
%     [blobRemove,CCfilt] = blobFilt(blobinfo,CC,BW)
% Inputs:
%     blobinfo,         Table with information of detected objects
%     CC,               Connected components, returned as a structure
%     BW,               Binary image of the segmented objects

% Outputs:
%     blobRemove,       Filtered table with information of detected objects
%% Remove top and bottom centroids
% Remove objects that were at the first 3 and last 3 images
stackHeight= size(BW,3);
% Remove centroids around the edge centroids from plane 3.5 and down
remv1 = blobinfo.Centroid(:,3)<=3;
% Detect the length of pictures and remove the last 3
remv2 = blobinfo.Centroid(:,3)>(stackHeight-3);
% Removing objects
remv =logical(remv1+remv2);
CC.PixelIdxList(remv)=[] ;
blobinfo(remv,:) = [];
blobRemove = blobinfo;
CCfilt = CC;
end


function filtVoxels =filtEstimateFnc(radius,xyRes,z)
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

% Equation of a sphere
SphereVol = 4/3*pi*radius^3;
voxel = xyRes/1000*xyRes/1000*z/1000; % voxels required for a perfect sphere
% estimate number of necessary voxels
filtVoxels = round(SphereVol./voxel);

%%%%%%%%%%%%%%%% Example %%%%%%%%%%%%%%%%%%%%%
% SphereVol = 4/3*pi*r^3
% Estimate volume with a radius of 3µm
% SphereVol = 4/3*pi*3^3 = 113µm^3
% voxel of pixels = 0.272*0.272*0.8µm = 0.06µm^3/voxel
% Estimere antal voxel som er nødvendige med en radius af 3µm
% voxel = 113µm^3./0.06µm^3/voxel = 1884
% tr= (4/3*pi*7^3)/0.06;
% length(find(blobfilt.Volume>tr))
end

function [folder,folderData] =folderGenerateFnc(folder1,folder2,folder3)
% Create the save folder of images
% Sintax:
%     [folder,folderData] =folderGenerateFnc(folder1,folder2,folder3,fileName)
% Inputs:
%     folder1,          First folder name
%     folder2,          Second folder name
%     folder3,          Third folder name
%     fileName,         Name of predticted images from UNetDense

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
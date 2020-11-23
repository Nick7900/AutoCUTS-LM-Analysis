% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% Analysing pyramidal cells
VISUALIZE = false;
SMOOTH =false; % smooth the 3D-reconstruction image
CROP = false;
SAVEFILES = true;
res = 272   ; % output resolution of data
xy = 272; % xy resolution (nm)
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
runProgram(folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,name,side,VISUALIZE,SMOOTH,CROP,SAVEFILES);

%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
function  runProgram(folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,name,side,VISUALIZE,SMOOTH,CROP,SAVEFILES)
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
%     CROP,        Crop area again if necessary
%     SAVEFILES,   Save data

%%%%%%%%%%%%%%%%% Start program
radius = []; % filtering noise objects away.
[BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,~] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius);
[blobfilt,CCfilt] =blobFnc(BW,z,res);
%%
% Start k-means to filter out small neurons
disp('Kmeans')
distThresh = 1.5; % ratio of the kmeans
[idx,C,kData,roundNeurons,pyramidalNeurons] =kmeans2Fnc(blobfilt,folder,distThresh);
[roundNeurons2] =kmeans3Fnc(blobfilt,roundNeurons,xy,res,folder);
roundNeurons =roundNeurons2;
kMeanOut = struct('idx',idx,'C',C,'kData',kData,'distThresh',distThresh,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
clear idx kData distThresh
%% Visualize and save files
if VISUALIZE == 0
    [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side);
    if SAVEFILES
        %%%%%%%%%%%%%%%%%%% Save images of binary segmentation
        for i = 1:size(pyramideNeuronImage,3)
            disp(['Save Image number ',num2str(i)])
            imwrite(pyramideNeuronImage(:,:,i),[folderImageSave,'/',fileNameSave,num2str(i),'.png'],'png');
            imwrite(roundNeuronImage(:,:,i),[folderImageSave2,'/',fileNameSave2,num2str(i),'.png'],'png');
        end
        disp('All images have been saved')
    end
else

    disp('Visualize 3D-reconstruction of pyramidal cells')
    [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side,SAVEFILES,SMOOTH);
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

function [blobfilt,CCfilt] =blobFnc(BW,z,res)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobfilt,CCfilt] =blobFnc(BW,z,res)
% Inputs:
%     BW,               Binary image stack
%     z,                z pixel resolution of images
%     res,              Output voxel resolution

% Outputs:
%     blobfilt,   Table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure

disp('Detect connections')
CC = bwconncomp(BW);
blobinfo = regionprops3(CC,'BoundingBox','Centroid','Volume','Orientation','VoxelList','PrincipalAxisLength','SurfaceArea');
[blobRemove, CCfilt]= blobFilt(blobinfo,CC,BW);
[blobfilt,CCfilt] = blobParameters(blobRemove,CCfilt,z,res);
end

function [BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,cRect] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius)
% Get path to read files
% Sintax:
%     [BW,folderImageSave,fileNameSave,folderImageSave2,fileNameSave2,folder,cRect] =readImages(CROP,folder1,folder2,folder3,fileName,ImgType,res,xy,z,buffer,radius)
% Inputs:
%     folder1,      First folder name
%     folder2,      Second folder name
%     folder3,      Third folder name
%     fileName,     Name of image that we read from
%     ImgType,      Image type (eg. TIF, png, JPEG)
%     res,          Output voxel resolution
%     xy,           xy pixel resolution of images
%     z,            z resolution (distance between each image)
%     buffer,       starting point from analysis the images (eg. buffer 20
%                   => analysise from image 20)
%     radius,       Artefact size based on a defined radius of a perfect sphere.

% Outputs:
%     srcFiles,     Structure with file informations
%     path,         path directory of where images should be read from

%%%%%%%%%%%%%%%%%% Sitting variables

%%%%%%%%%%%% filter small cells away
% Filter noise less than 2.1904 µm^3 => 37 voxels default settings
filtVoxels =filtEstimateFnc(radius,res,z);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create save image folder %%%%%%%%%%

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


%% Read Image
[I,Ic,cRect] = readImgFnc(folder1,folder2,folder3,fileName,ImgType,res,xy,buffer,CROP);

%% Filter and fill the holes of the images
if ~isempty(Ic) == 1
    disp('Filter Crop image')
    [BW]= filtFillFnc(Ic, filtVoxels);
else
    disp('Filter Original image')
    [BW]= filtFillFnc(I, filtVoxels);
end
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
% Default radius is 0.8 µm
if isempty(radius)
    radius =0.8; % decide a radius in µm
    SphereVol = 4/3*pi*radius^3;
    voxel = res/1000*res/1000*z/1000;
    filtVoxels = round(SphereVol./voxel);
else
    SphereVol = 4/3*pi*radius^3;
    voxel = res/1000*res/1000*z/1000;
    % estimate number of necessary voxels
    filtVoxels = round(SphereVol./voxel);
end

% SphereVol = 4/3*pi*r^3
% Estimate volume with a radius of 3µm
% SphereVol = 4/3*pi*3^3 = 113µm^3
% voxel of pixels = 0.272*0.272*0.8µm = 0.06µm^3/voxel
% Estimere antal voxel som er nødvendige med en radius af 3µm
% voxel = 113µm^3./0.06µm^3/voxel = 1884
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

Img =imread([path,fileName,num2str(1+buffer),'.',ImgType]);
[~, ~, numberOfColorChannels] = size(Img);

% Change to desired image resolution
if res==xy
    if numberOfColorChannels > 1
        % Convert to grayscale
        ImgGray =rgb2gray(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),N);
        I(:,:,1) = imbinarize(ImgGray);
        %%
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N)])
            I_ori = (rgb2gray((imread([path,fileName,num2str(i),'.',ImgType]))));
            I(:,:,i) =imbinarize(I_ori);
        end
        clc
        disp([num2str(N),' binary images are loaded '])
    else
        ImgGray =(Img);
        %%%%%%%%%%%%%%% Create the image matrix for volume data
        I =false(size(ImgGray,1),size(ImgGray,2),N);
        I(:,:,1) = (ImgGray);
        
        for i = 2+buffer:N+buffer
            disp(['Load image number ',num2str(i),' out of ',num2str(N+buffer)])
            I_ori= (((imread([path,fileName,num2str(i),'.',ImgType]))));
            I(:,:,i-buffer) =(I_ori);
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

function [idx,C,kData,roundNeurons,pyramidalNeurons] =kmeans2Fnc(blobfilt,folder,distThresh)
% Seperating the dataset between estimated radius and farthest voxel
% from centroid by using kmeans
% Sintax:
%     [idx,C,kData,roundNeurons,pyramidalNeurons] =kmeans2Fnc(blobfilt,folder,distThresh)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     folder,           folder name
%     distThresh,       Defined threshold ratio between the farthest voxel
%                       from centroid and estimated radius

% Outputs:
%     idx,              Cluster indices, returned as a numeric column vector.
%     C,                Cluster centroid locations, returned as a numeric matrix
%     kData,            Data used for kmeans
%     roundNeurons,     Detected round neurons
%     pyramidalNeurons, Detected pyramidal cells
VISUALIZE =1;

% Define variables
clus = 2;
blobfiltDisk2 = blobfilt(blobfilt.DistRatio<distThresh, :);
blobfiltDisk = blobfiltDisk2(blobfiltDisk2.Radius<5, :);
% K means%%%%%%%%%%%%%%%%%%%%%%%%%%%
x= [blobfiltDisk.Radius];
y= [blobfiltDisk.DistMax];
opts = statset('Display','final');
kData = [x y];
[idx,C] = kmeans(kData,clus,'dist','sqeuclidean',...
    'Replicates',10,'Options',opts);
% 'sqeuclidean' (default) | 'cityblock' | 'cosine' | 'correlation' | 'hamming'

roundC =C(:,2);
group =find(roundC==min(C(:,2)));
if VISUALIZE
    %# show points and clusters (color-coded)
    
    
    if size(C,1)==2
        figure, hold on
        if group==1
            plot(kData(idx==1,1),kData(idx==1,2),'r.','MarkerSize',10)
            plot(kData(idx==2,1),kData(idx==2,2),'b.','MarkerSize',10)
        else
            plot(kData(idx==1,1),kData(idx==1,2),'b.','MarkerSize',10)
            plot(kData(idx==2,1),kData(idx==2,2),'r.','MarkerSize',10)
        end
        plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',2)
        xlabel('Radius'), ylabel('DistMax')
        legend('Cluster 1','Cluster 2','Centroids',...
            'Location','NE')
        title 'Cluster Assignments and Centroids'
        hold off
    elseif size(C,1)==3
        figure, hold on
        plot(kData(idx==1,1),kData(idx==1,2),'r.','MarkerSize',10)
        plot(kData(idx==2,1),kData(idx==2,2),'b.','MarkerSize',10)
        plot(kData(idx==3,1),kData(idx==3,2),'g.','MarkerSize',10)
        plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',2)
        xlabel('Volume'), ylabel('Sphericity')
        legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
            'Location','NE')
        title 'Cluster Assignments and Centroids'
        hold off
        
    elseif size(C,1)==4
        figure, hold on
        plot(kData(idx==1,1),kData(idx==1,2),'r.','MarkerSize',10)
        plot(kData(idx==2,1),kData(idx==2,2),'b.','MarkerSize',10)
        plot(kData(idx==3,1),kData(idx==3,2),'g.','MarkerSize',10)
        plot(kData(idx==4,1),kData(idx==4,2),'y.','MarkerSize',10)
        plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',2)
        hold off
        xlabel('Volume'), ylabel('Sphericity')
        legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
            'Location','NE')
        title 'Cluster Assignments and Centroids'
        
    else
        clr = lines(clus);
        figure, hold on
        scatter(kData(:,1), kData(:,2), 36, clr(idx,:), 'Marker','.')
        scatter(C(:,1), C(:,2), 100, clr, 'Marker','o', 'LineWidth',3)
        hold off
        xlabel('Volume'), ylabel('Sphericity')
        title ('Cluster Assignments and Centroids')
    end
    
    
end
saveas(gcf,fullfile(folder, 'kMeansplot'),'m')
saveas(gcf,fullfile(folder, 'kMeansplot'),'png')
hold off
maxSpher = zeros(clus,1);
%     Try to find the maximum volume
for i = 1:clus
    maxSpher(i)= max(y(idx==i,:));
end
spC =C(:,2);
roundNeurons =idx==find(spC==min(C(:,2))); % Detect the cluster with the highest specity

%
% Detect minimum radius for pyramidal cells
blobfiltDisk3 = blobfiltDisk;
blobfiltDisk3((roundNeurons),:) = []; % Radius of cells that are not pyramidal cells
minRadius = min(blobfiltDisk3.Radius);

% Add remaining logical numbers
roundNeurons = logical([roundNeurons;zeros(size(blobfilt,1)-size(roundNeurons,1),1)]);
%%%%%%%%%%%%%%% Filter small volumes away as well
remvSmall =(blobfilt.BoundingBox(:,6)<=3);
remvSmall2 =(blobfilt.Radius<minRadius);
roundNeurons =logical(roundNeurons+remvSmall+remvSmall2);
pyramidalNeurons =~logical(roundNeurons); % the rest are round cells
end

function [roundNeurons] =kmeans3Fnc(blobfilt,roundNeurons,xy,res,folder)
% Seperating the dataset between estimated boundingbox size and estimated sphericity
% Sintax:
%     [roundNeurons] =kmeans3New3Fnc(blobfilt,roundNeurons,xy,res,folder)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     roundNeurons,     Detected round neurons
%     xy,               xy pixel resolution of images
%     res,              Output voxel resolution
%     folder,           folder name

% Outputs:
%     roundNeurons,     Detected round neurons

VISUALIZE =1;
cluster = 2; % number of clusters
% Wont estimate kmean3 if the stack doesnt include more than 100 images
stack_z = max(blobfilt.BoundingBox(:,3));  
limit = round(110*(xy/res)); %length of 110*0.272 nm % 30 µm                 
limit1 = round(147*(xy/res)); %length of 147*0.272 nm % 40 µm             

blobfiltDisk2 = blobfilt(blobfilt.BoundingBox(:,4)>limit1, :);
blobfiltDisk = blobfiltDisk2(blobfiltDisk2.BoundingBox(:,5)>limit, :);

idxBound = find(blobfilt.BoundingBox(:,4)>limit1&blobfilt.BoundingBox(:,5)>limit);
if size(blobfiltDisk,1)<3
    blobfiltDisk = blobfiltDisk2;
    idxBound = find(blobfilt.BoundingBox(:,4)>limit1);
end

if size(blobfiltDisk,1)<3||stack_z<100
     disp('<strong>No big random blobs</strong>')
else
    % K means%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xVal= blobfiltDisk.BoundingBox(:,4);
yVal= blobfiltDisk.BoundingBox(:,5);
    zVal= blobfiltDisk.Sphericity;
    [normX] =normFnc(xVal);
    [normY] =normFnc(yVal);
    [normZ] =normFnc(zVal);
    X = [normX normY normZ];
    
    opts = statset('Display','final');
    [idx,C] = kmeans(X,cluster,'Distance','cityblock',...
        'Replicates',100,'Options',opts);
    spC =C(:,3);
    %     abs(diff(C(:,3)))
    
    if sum(blobfiltDisk.Sphericity<0.23)>0 && abs(diff(C(:,3)))>0.25
        group =find(spC==min(C(:,3)));
        if VISUALIZE
            %# show points and clusters (color-coded)
            clr = lines(cluster);
            hFig = figure();
            axh = axes('Parent', hFig);
            hold(axh, 'all');
            if group ==1
                h1 =scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),50,'b.','Marker','o', 'LineWidth',2);
                h2 =scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),50,'r.','Marker','o', 'LineWidth',2);
            else
                h1 =scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),50,'b.','Marker','o', 'LineWidth',2);
                h2 =scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),50,'r.','Marker','o', 'LineWidth',2);
            end
            h3 =scatter3(C(:,1), C(:,2), C(:,3), 100, clr,'k.', 'Marker','x', 'LineWidth',3);
            hold off
            view(42,16), axis vis3d, box on, rotate3d on
            xlabel('X'), ylabel('Y'), zlabel('Z')
            title(['Random blobs in Cluster: ',num2str(group)])
            xticks((linspace(min(X(:,1)),max(X(:,1)),5)));
            yticks((linspace(min(X(:,2)),max(X(:,2)),5)));
            zticks((linspace(min(X(:,3)),max(X(:,3)),4)));
            legend(axh, [h1,h2,h3], {'Cluster 1', 'Cluster 2','Centroids'});
        end
        maxBox = zeros(cluster,1);
        %     Try to find the maximum volume
        for i = 1:cluster
            maxBox(i)= max(zVal(idx==i,:));
        end
        
        saveas(gcf,fullfile(folder, 'kMeansBoundplot'),'m')
        saveas(gcf,fullfile(folder, 'kMeansBoundplot'),'png')
        roundNeuronsBound =idx==find(spC==min(C(:,3))); % Detect the cluster with the lowest specity
        roundNeurons(idxBound(roundNeuronsBound)) =1;
        
    else
        disp('<strong>No big random blobs</strong>')
    end
end
%%
end

function [norm_data] =normFnc(x)
% Normalize data
% Sintax:
%     [norm_data] =normFnc(x)
% Inputs:
%     x,              Input data
% Outputs:
%     norm_data,      Normalized data
%% Normalize data
minVal = min(x);
maxVal = max(x);
norm_data = (x - minVal) / ( maxVal - minVal); %between 0 and 1
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
remv = blobinfo.Centroid(:,3)<3.5;
remv1 =blobinfo.BoundingBox(:,6)>35;
remv =logical(remv+remv1);
CC.PixelIdxList(remv)=[] ;
blobinfo(remv,:) = [];
% Detect the length of pictures and remove the last 3
lengthBlob= size(BW,3);
remv = blobinfo.Centroid(:,3)>(lengthBlob-3.5);
blobinfo(remv,:) = [];
CC.PixelIdxList(remv)=[] ;

% find(blobfilt.BoundingBox(:,6)<3);

%% Find objects where the orientation is 0 in the z-axis which mean it is only a single plane object
idxOri = find(blobinfo.Orientation(:,3)==0);
meanVol =mean(blobinfo.Volume(:));
medianVol =median(blobinfo.Volume(:));
% find the values of specific index and convert them to 0
if meanVol<medianVol
    biOri =(blobinfo.Volume(idxOri)<meanVol);
    removeBlob=biOri.*idxOri;
    % Check if there are some who are empty or not
    if find(removeBlob==0)>1
        removeBlob(removeBlob==0)=[];
    end
    blobinfo(removeBlob,:)=[];
    CC.PixelIdxList(:,removeBlob')=[];
    CC.NumObjects = size(CC.PixelIdxList,2);
else
    % find the values of specific index
    biOri =(blobinfo.Volume(idxOri)<medianVol);
    removeBlob=biOri.*idxOri;
    % Check if there are some who are empty or not and remove them
    %     if find(removeBlob==0)>1
    %         removeBlob(removeBlob==0)=[];
    %     end
    removeBlob(removeBlob==0)=[];
    blobinfo(removeBlob,:)=[];
    %     blobinfo(logical(removeBlob),:)=[];
    CC.PixelIdxList(removeBlob)=[] ;
    %     CC.PixelIdxList(logical(removeBlob))=[] ;
    CC.NumObjects = size(CC.PixelIdxList,2);
end

blobRemove = blobinfo;
CCfilt = CC;
end

function [blobfiltout,CCfilt] = blobParameters(blobRemove,CCfilt,z,res)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobfiltout,CCfilt] = blobParameters(blobRemove,CCfilt,z,res)
% Inputs:
%     blobRemove,       Filtered table with information of detected objects
%     CCfilt,           Filtered connected components, returned as a structure
%     z,                z resolution (distance between each image)
%     res,              Output voxel resolution

% Outputs:
%     blobfiltout,      Table with information of detected objects with
%                       added parameters
%     CCfilt,           Filtered Connected components, returned as a structure

%% Estimate the longest distance between centroid and all voxels
centroids =blobRemove.Centroid;
distMax =zeros(size(blobRemove,1),1);
diskLine =zeros(size(blobRemove,1),3);

for i = 1:size(blobRemove,1)
    distArray = zeros(size(blobRemove.VoxelList{i,1},1),1); %distance between centroid and voxels
    centroid1=centroids(i,:);
    voxelList =blobRemove.VoxelList{i,1};
    for j = 1:size(blobRemove.VoxelList{i,1},1)
        distance_between_points = sqrt(sum((centroid1 - (voxelList(j,:))).^2)); % distance
        distArray(j) = distance_between_points;
    end
    % Find the farthest voxel from centroid
    distMax(i)= max(distArray);
    idxMax = distArray==distMax(i); % locical array of max distance
    distOne = find(idxMax, 1, 'first'); % if there more than one index with the length, then use the first one
    diskLine(i,:)=voxelList(distOne,:);
    
end
%% % Estimate the distance between centroid and all voxels
% distance_between_points = sqrt(((centroid1(1) - voxelList(:,1)).^2 + ((centroid1(2) - voxelList(:,2)).^2)+((centroid1(3) - voxelList(:,3)).^2)));
blobRemove.DistMax = distMax*(res/1000);% Convert to micrometer
blobRemove.diskLine = diskLine;
adjustAxisLength =blobRemove.PrincipalAxisLength;
adjustAxisLength(:,3)=adjustAxisLength(:,3).*z/res; % Convert z to xy resolution pr pixel
%https://math.wikia.org/wiki/Ellipsoidal_quadratic_mean_radius
% Convert the length of axis to µm
adjustMicrometer =adjustAxisLength* (res/1000); % Convert length to micrometer
blobRemove.PrincipalAxisLengthAdjust = adjustMicrometer;
diameters = sqrt(sum(blobRemove.PrincipalAxisLengthAdjust.^2,2)/3); % average/mean length of principalAxis
blobRemove.Radius =diameters./2; % Convert the average length to radius
blobRemove.DistRatio =blobRemove.DistMax./blobRemove.Radius;

%% sort rows && Output
[blobfiltout, idxSort] = sortrows(blobRemove,size(blobRemove,2),'ascend');
%%%% Indext CCfilt
CCfilt.PixelIdxList=CCfilt.PixelIdxList(idxSort);
V=blobfiltout.Volume;
blobfiltout.Sphericity = (pi^(1/3)*(6*V).^(2/3))./blobfiltout.SurfaceArea;%Sphericity
%%%%% Estimage surface of ellipse %%%%%%%%%%%%% Not used
% ellipsMicro =blobfiltout.PrincipalAxisLengthAdjust;
% total length has to be divided in 2 to get the radius of all 3 axis
% p =1.6075;
% a = ellipsMicro(:,1)/2;
% b = ellipsMicro(:,2)/2;
% c = ellipsMicro(:,3)/2;
% SA = 4*pi*(((a.*b).^p+(a.*c).^p+(b.*c).^p)/(3)).^(1/p);
% volEllipse = 4/3.*pi.*(a.*2).*(b.*2).*(c.*2)
end

function [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side,SAVEFILES,SMOOTH)
% Save data of the 3D-reconstruction of pyramidal cells
% Sintax:
%     [pyramideNeuronImage,roundNeuronImage] =visSaveDataReconstructionFnc(blobfilt,CCfilt,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side,SAVEFILES,SMOOTH)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure
%     BW,               Binary image stack
%     pyramidalNeurons, Detected pyramidal neurons
%     roundNeurons,     Detected round neurons + artefacts
%     res,              Output voxel resolution
%     xy,               xy pixel resolution of images
%     z,                z pixel resolution of images
%     folder,           folder name
%     kMeanOut,         kMean output data
%     name,             name of subject/folder that will be analysed
%     side,             Pial surface direction: left side 0, right side 1
%     SAVEFILES,        Save data
%     SMOOTH,           Smooth 3D-reconstruction images

% Outputs:
%     pyramideNeuronImage,  Image stack of Pyramidal cells
%     roundNeuronImage,     Image stack of round neurons + artefacts

%% Define sampling window
window = [size(BW,2), size(BW,1),round(size(BW,3).*z/res)];  %[xLeft, yTop, length, width, hieght].
z_heght = size(BW,3);
clear BW;

%% Visualize neurons
roundNeurons2 =logical(roundNeurons);
pyramidalNeurons2 = logical(pyramidalNeurons.*(~roundNeurons2));
neuronSmall = (find(roundNeurons2));
neuronBig =(find(pyramidalNeurons2));
blobfiltPyramid = blobfilt;
blobfiltPyramid(~(pyramidalNeurons2),:) = [];

%%

if SMOOTH
    labeledImage = labelmatrix(CCfilt);
    d = datestr(datetime('today'));
    % Only find the big and small volumes
    smallNeuronImage = ismember(labeledImage, neuronSmall);
    pyramideNeuronImage = ismember(labeledImage, neuronBig);
    
    % centroids
    cenAdRound = blobfilt.Centroid(neuronSmall,:); % finding the centriods for small neurons
    cenAdPyra = blobfilt.Centroid(neuronBig,:); % finding the centriods for pyramid neurons
    diskLineHigh = blobfilt.diskLine(neuronBig,:);
    
    %%%%%%%%%%%% smooth the surface of the neurons for visualization
    dataRound = smooth3(smallNeuronImage,'box',3);
    dataPyra = smooth3(pyramideNeuronImage,'box',3);
    
    %%%%%% Visualize the 3D-reconstruction
    figure;
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
    xticks(linspace(1,size(labeledImage,2),8));
    xticklabels(round(linspace(1,window(1),8)*res/1000));
    yticks(linspace(1,size(labeledImage,1),8));
    yticklabels(round(linspace(1,window(2),8)*res/1000));
    zticks(linspace(1,size(labeledImage,3),3));
    zticklabels(round(linspace(1,window(3),3)*res/1000));
    %         set(gca,'DataAspectRatio', [1 1 1]);
    title(['Pyramid cells crop sequence ',d],'FontSize',14);
    saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCells'),'png')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 3D orientation plot
    % Estimate orientation of each pyramidal cell
    blobfiltPyramid.theta =oriEstimateFnc(window,blobfiltPyramid,side);
    % Mark orientation above 90 degrees
    filterOri = blobfiltPyramid.theta>90;
    
    % Visualize orientation of pyramidal cells
    for i =1:size(cenAdPyra,1)
        P1 = cenAdPyra(i,:);
        P2 = diskLineHigh(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','yellow', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    % Visualize orientations above 90 degrees with red lines
    diskLineOriFilt = blobfiltPyramid.diskLine(filterOri,:);
    cenOriFilt = blobfiltPyramid.Centroid(filterOri,:);
    hold on
    for i =1:size(cenOriFilt,1)
        P1 = cenOriFilt(i,:);
        P2 = diskLineOriFilt(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','red', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri'),'png')
    close all
else
    labeledImage = labelmatrix(CCfilt);
    d = datestr(datetime('today'));
    % Only find the big and small volumes
    dataRound = ismember(labeledImage, neuronSmall);
    dataPyra = ismember(labeledImage, neuronBig);
    
    % centroids
    cenAdRound = blobfilt.Centroid(neuronSmall,:); % finding the centriods for small neurons
    cenAdPyra = blobfilt.Centroid(neuronBig,:); % finding the centriods for pyramid neurons
    diskLineHigh = blobfilt.diskLine(neuronBig,:);
    
    %%%%%% Visualize the 3D-reconstruction
    figure;
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
    xticks(linspace(1,size(labeledImage,2),8));
    xticklabels(round(linspace(1,window(1),8)*res/1000));
    yticks(linspace(1,size(labeledImage,1),8));
    yticklabels(round(linspace(1,window(2),8)*res/1000));
    zticks(linspace(1,size(labeledImage,3),3));
    zticklabels(round(linspace(1,window(3),3)*res/1000));
    %         set(gca,'DataAspectRatio', [1 1 1]);
    title(['Pyramid cells crop sequence ',d],'FontSize',14);
    saveas(gcf,fullfile(folder, '3D_pyramidalCells_raw'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCells_raw'),'png')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 3D orientation plot
    % Estimate orientation of each pyramidal cell
    blobfiltPyramid.theta =oriEstimateFnc(window,blobfiltPyramid,side);
    % Mark orientation above 90 degrees
    filterOri = blobfiltPyramid.theta>90;
    
    % Visualize orientation of pyramidal cells
    for i =1:size(cenAdPyra,1)
        P1 = cenAdPyra(i,:);
        P2 = diskLineHigh(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','yellow', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    % Visualize orientations above 90 degrees with red lines
    diskLineOriFilt = blobfiltPyramid.diskLine(filterOri,:);
    cenOriFilt = blobfiltPyramid.Centroid(filterOri,:);
    hold on
    for i =1:size(cenOriFilt,1)
        P1 = cenOriFilt(i,:);
        P2 = diskLineOriFilt(i,:);
        pts = [P1; P2];
        line(pts(:,1), pts(:,2), pts(:,3),'Color','red', 'LineWidth',3)
        plot3(pts(:,1), pts(:,2), pts(:,3),'k')
        hold on
    end
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri_raw'),'m')
    saveas(gcf,fullfile(folder, '3D_pyramidalCellsOri_raw'),'png')
    close all
end
if SAVEFILES==1
    %% Estimate the two centroids for round and pyramidal cells
    cenLow = blobfilt.Centroid(roundNeurons2,:);
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
    disp('Save variables')
    % kMeanOut = struct('Index',idx,'Clusters',C,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
    save(fullfile(folder, 'kMeanOut.mat'), 'kMeanOut')
    save(fullfile(folder, 'window.mat'), 'window')
    save(fullfile(folder, 'blobfilt.mat'), 'blobfilt', '-v7.3')
    save(fullfile(folder, 'cenAdLow.mat'), 'cenAdLow')
    save(fullfile(folder, 'cenAdHigh.mat'), 'cenAdHigh')
    %%%%%%%%%% Save segmentation output
    save(fullfile(folder, 'blobfiltPyramid.mat'), 'blobfiltPyramid', '-v7.3')
    my_directory =fullfile(folder, ['centroidPyramidal_',name,'.xlsx']);
    xlswrite(my_directory,cenAdHigh); % Excel file of coordinates
    clear CCfilt kMeanOut blobfilt cenAdLow cenAdHigh blobfiltPyramid roundNeurons roundNeurons2
    %% Recontruct 3D image from labels bwconnected image
    disp('Reconstruct Pyramidal cell images')
    ImgSize =window(2)*window(1)*z_heght;
    newImgPyra = false(ImgSize,1);
    for i = 1:size(CCfiltAdjust.PixelIdxList,2)
        sizeCC =size(CCfiltAdjust.PixelIdxList{1, i},1);
        for j = 1:sizeCC
            newImgPyra(CCfiltAdjust.PixelIdxList{1, i}(j)) =CCfiltAdjust.PixelIdxList{1, i}(j);
        end
    end
    pyramideNeuronImage = reshape(newImgPyra,[window(2),window(1),z_heght]);
    save(fullfile(folder, 'pyramideNeuronImage.mat'), 'pyramideNeuronImage', '-v7.3')
    
    disp('Reconstruct Round cell images')
    ImgSize =window(2)*window(1)*z_heght;
    newImgRound = false(ImgSize,1);
    for i = 1:size(CCfiltRound.PixelIdxList,2)
        sizeCC =size(CCfiltRound.PixelIdxList{1, i},1);
        for j = 1:sizeCC
            newImgRound(CCfiltRound.PixelIdxList{1, i}(j)) =CCfiltRound.PixelIdxList{1, i}(j);
        end
    end
    roundNeuronImage = reshape(newImgRound,[window(2),window(1),z_heght]);
    save(fullfile(folder, 'roundNeuronImage.mat'), 'roundNeuronImage', '-v7.3')
end

end


function [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,res,xy,z,folder,kMeanOut,name,side)   
% Save data of the 3D-reconstruction of pyramidal cells
% Sintax:
%     [pyramideNeuronImage,roundNeuronImage] =saveDataReconstructionFnc(blobfilt,CCfilt,BW,pyramidalNeurons,roundNeurons,xy,z,folder,kMeanOut,name,side)
% Inputs:
%     blobfilt,         Filtered table with information of detected objects
%     CCfilt,           Filtered Connected components, returned as a structure
%     BW,               Binary image stack
%     pyramidalNeurons, Detected pyramidal neurons
%     roundNeurons,     Detected round neurons + artefacts
%     res,              Output voxel resolution
%     xy,               xy pixel resolution of images
%     z,                z pixel resolution of images
%     folder,           folder name
%     kMeanOut,         kMean output data
%     name,             name of subject/folder that will be analysed
%     side,             Pial surface direction: left side 0, right side 1


% Outputs:
%     pyramideNeuronImage,  Image stack of Pyramidal cells
%     roundNeuronImage,     Image stack of round neurons + artefacts
%% Define sampling window
window = [size(BW,2), size(BW,1),round(size(BW,3).*z/res)];  %[xLeft, yTop, length, width, hieght].
z_heght = size(BW,3);
clear BW;
%% Remove small neurons
roundNeurons2 =logical(roundNeurons);
pyramidalNeurons = logical(pyramidalNeurons.*(~roundNeurons));
pyramidalNeurons = logical(pyramidalNeurons.*(~roundNeurons2));
blobfiltPyramid = blobfilt;
blobfiltPyramid(~(pyramidalNeurons),:) = [];
% Estimate orientation of each pyramidal cell
blobfiltPyramid.theta =oriEstimateFnc(window,blobfiltPyramid,side);
%% Estimate the two centroids for round and pyramidal cells
cenLow = blobfilt.Centroid(roundNeurons2,:);
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
disp('Save variables')
% kMeanOut = struct('Index',idx,'Clusters',C,'RoundNeurons',roundNeurons,'PyramidalNeurons',pyramidalNeurons);
save(fullfile(folder, 'kMeanOut.mat'), 'kMeanOut')
save(fullfile(folder, 'window.mat'), 'window')
save(fullfile(folder, 'blobfilt.mat'), 'blobfilt', '-v7.3')
save(fullfile(folder, 'cenAdLow.mat'), 'cenAdLow')
save(fullfile(folder, 'cenAdHigh.mat'), 'cenAdHigh')
%%%%%%%%%% Save segmentation output
save(fullfile(folder, 'blobfiltPyramid.mat'), 'blobfiltPyramid', '-v7.3')
my_directory =fullfile(folder, ['centroidPyramidal_',name,'.xlsx']);
xlswrite(my_directory,cenAdHigh); % Out of memory when you write
clear CCfilt kMeanOut blobfilt cenAdLow cenAdHigh blobfiltPyramid roundNeurons roundNeurons2
%% Recontruct 3D image from labels bwconnected image
disp('Reconstruct pyramide images')
ImgSize =window(2)*window(1)*z_heght;
newImgPyra = false(ImgSize,1);
for i = 1:size(CCfiltAdjust.PixelIdxList,2)
    sizeCC =size(CCfiltAdjust.PixelIdxList{1, i},1);
    for j = 1:sizeCC
        newImgPyra(CCfiltAdjust.PixelIdxList{1, i}(j)) =CCfiltAdjust.PixelIdxList{1, i}(j);
    end
end
pyramideNeuronImage = reshape(newImgPyra,[window(2),window(1),z_heght]);
save(fullfile(folder, 'pyramideNeuronImage.mat'), 'pyramideNeuronImage', '-v7.3')

disp('Reconstruct Round images')
ImgSize =window(2)*window(1)*z_heght;
newImgRound = false(ImgSize,1);
for i = 1:size(CCfiltRound.PixelIdxList,2)
    sizeCC =size(CCfiltRound.PixelIdxList{1, i},1);
    for j = 1:sizeCC
        newImgRound(CCfiltRound.PixelIdxList{1, i}(j)) =CCfiltRound.PixelIdxList{1, i}(j);
    end
end
roundNeuronImage = reshape(newImgRound,[window(2),window(1),z_heght]);
save(fullfile(folder, 'roundNeuronImage.mat'), 'roundNeuronImage', '-v7.3')
end

function theta =oriEstimateFnc(window,blobfiltPyramid,side)
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
    c = blobfiltPyramid.Centroid; % position of centroid
    d =blobfiltPyramid.diskLine; % positio of voxel farthest away from centroid
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
    c = blobfiltPyramid.Centroid; % position of centroid
    d =blobfiltPyramid.diskLine; % position of voxel farthest away from centroid
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

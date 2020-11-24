% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% Objectwise Validation
folder1 ='Manuel_prepared';
folder2 ='Prediction_prepared';
folder3 ='Original_prepared';
fileName = 'ori';
ImgType = 'png';
%% Start program
% Read Image
[I_mark] = readImg(folder1,fileName,ImgType);
[I_pred] = readImg(folder2,fileName,ImgType);
[I_ori] = readOriImg(folder3,fileName,ImgType);
disp('Analyse sections')
%% Filter and fill the holes of the images
% Find connected components
[blobMark,idxCenMark] =blobFnc(I_mark);
[blobPred,idxCenPred] =blobFnc(I_pred);
%%%%%%% Generate overlay folder
[foldDataMark] =createOverlayFolderFnc(folder1);
[foldDataPred,foldDataFP] =createOverlayFolderFnc(folder2);

%% Save variables and calculation
saveVariablesFnc(I_mark,I_pred,I_ori,idxCenMark,idxCenPred,blobMark,blobPred,foldDataPred,foldDataFP)
%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
% Create the varibles for blob output
function [I_mark] = readImg(folder1,fileName,ImgType)
% Reading image stack in specified folder
% Sintax:
%     [I_mark] = readImg(folder1,fileName,ImgType)
% Inputs:
%     folder1,     First folder name
%     fileName,    Name of image file
%     ImgType,     Image type (eg. TIF, png, JPEG)

% Outputs:
%     I_mark,      Image output

srcFiles = dir(fullfile([folder1,'/*.',ImgType]));
path=srcFiles.folder;
path=[path '\'];

N =size(srcFiles,1);
% read the first image
I = imread([path,fileName,' (',num2str(1),').',ImgType]);
%% Convert image to grayimage and load images in folder
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[~, ~, numberOfColorChannels] = size(I);
if numberOfColorChannels > 1
    % It's not really gray scale like we expected - it's color.
    I_mark = (rgb2gray(I));
    I_mark = repmat(I_mark,1,1,N);
    for i = 2:N
        disp(['Load image number ',num2str(i),' out of ',num2str(N)])
        I_mark(:,:,i) = (rgb2gray(imread([path,fileName,' (',num2str(i),').',ImgType])));
    end
    clc
    disp([num2str(N),' binary images are loaded '])
    
else
    I_mark = imbinarize(imread(I));
    I_mark = repmat(I_mark,1,1,N);
    for i = 2:N
        disp(['Load image number ',num2str(i),' out of ',num2str(N)])
        I_mark(:,:,i) = ((imread([path,fileName,' (',num2str(i),').',ImgType])));
    end
    
end

I_mark=imbinarize(I_mark);
clc
disp([num2str(N),' binary images are loaded '])
end


function [I_out] = readOriImg(folder3,fileName,ImgType)
% Reading image stack in specified folder
% Sintax:
%     [I_out] = readOriImg(folder3,fileName,ImgType)
% Inputs:
%     folder3,     folder name
%     fileName,    Name of image file
%     ImgType,     Image type (eg. TIF, png, JPEG)

% Outputs:
%     I_out,      Image output

srcFiles = dir(fullfile([folder3,'/*.',ImgType]));
path=srcFiles.folder;
path=[path '\'];

N =size(srcFiles,1);
% read the first image
I_out = imread([path,fileName,' (',num2str(1),').',ImgType]);

% I_out = imread([path,fileName',' (',num2str(1),').',ImgType]);
for i = 2:N
    disp(['Load image number ',num2str(i),' out of ',num2str(N)])
    I_out(:,:,i) = imread([path,fileName,' (',num2str(i),').',ImgType]);
    %     I_out(:,:,i) = uint8(((imread([path,fileName',' (',num2str(i),').',ImgType]))));
end
clc
disp([num2str(N),' original images are loaded ']);
end


function [blobInfo,idxCen] =blobFnc(I)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     [blobInfo,idxCen] =blobFnc(I)
% Inputs:
%     I,               Binary image stack

% Outputs:
%     blobInfo,         Table with information of detected objects
%     idxCen,           Index of the detected centrois 

% Detect connected components of the imagestack
CC = bwconncomp(I);
blobInfo = regionprops3(CC,'BoundingBox','Centroid');
% blobInfo = regionprops3(CC,'BoundingBox','Centroid','Volume','SurfaceArea');
blobInfo.CentroidRound = round(blobInfo.Centroid); % round centroids
blobInfo.Centroid =single(blobInfo.Centroid); % convert to single
blobInfo.BoundingBox =single(blobInfo.BoundingBox); % convert to single
% Detect connected compoents with at least three subsequent images   
remv1 =(blobInfo.BoundingBox(:,6)<=3); 
remv =logical(remv1);
blobInfo(remv,:) = []; % filter objects away
CC.PixelIdxList(remv)=[] ;
CC.NumObjects = size(CC.PixelIdxList,2);
%%
% Matrix of the centroids in each layer
% detect centroids
idxCen = false(size(blobInfo,1),size(I,3));
for i = 1:size(I,3)
    idxCen(:,i)=blobInfo.CentroidRound(:,3)==i;
end

% V=blobInfo.Volume;
% A=blobInfo.SurfaceArea;
% sp =(pi^(1/3)*(6*V).^(2/3))./A;%Sphericity
% blobInfo.Sphericity = single(sp);
end

function [foldDataPred,foldDataFP] =createOverlayFolderFnc(folder1)
% Create the save folder for the images
% Sintax:
%     [foldDataPred,foldDataFP] =createOverlayFolderFnc(folder2);
% Inputs:
%     folder1,     First folder name

% Outputs:
%     foldData,         Folder for predicted results
%     foldDataFP,       Folder for False positive results

srcFiles = dir(fullfile(folder1)); % Directory
[parentFolder] = fileparts(srcFiles(1).folder);

% Creating the folder overlay in the Result-folder
foldDataPred=('OverlayPrediction');
foldDataPred = sprintf(['%s/',foldDataPred,'%s'], parentFolder);
if ~exist(foldDataPred, 'dir')
    mkdir(foldDataPred);
end
foldDataFP=('OverlayFP');
foldDataFP = sprintf(['%s/',foldDataFP,'%s'], parentFolder);
if ~exist(foldDataFP, 'dir')
    mkdir(foldDataFP);
end
end

function saveVariablesFnc(I_mark,I_pred,I_ori,idxCenMark,idxCenPred,blobMark,blobPred,foldDataPred,foldDataFP)
% Measure properties of 3D volumetric objects from stacked images
% Sintax:
%     saveVariablesFnc(I_pred,idxCenMark,I_ori,blobMark,foldDataPred,I_mark,idxCenPred,blobPred,foldDataFP)
% Inputs:
%     I_mark,               Manuel mark images (Validationset)
%     I_pred,               UNetDense predicted output
%     I_ori,                Original grayscale images
%     idxCenMark,           Index of the manuel detected centrois 
%     idxCenPred,           Index of the predicted detected centrois 
%     blobMark,             Table with information of manuel detected objects
%     blobPred,             Table with information of predicted detected objects
%     foldDataPred,         Folder for predicted results
%     foldDataFP,           Folder for False positive results

%% Estimere TP and FN for predicted images from UNetDense
[tp,fn,~,~]  = sensitivityFnc(I,I_ori,idxCenMark,blobMark,foldDataPred);
[~,fp,~,~]  = precisionFnc(I_mark,I_ori,idxCenPred,blobPred,foldDataFP);
%%  Measures of the performance of a binary classification test of the model
TP=sum((tp));
FN=sum((fn)); % Neurons should be detected, but the program didnt detect any
FP=sum((fp)); % No neurons should be detected, but the program did detect a neuron
sensitivity = TP/(TP+FN);
precision = TP/(TP+FP);
F1 = 2/((sensitivity)^-1+(precision)^-1);

%%%%%%%%%%%%%%%%%% figure of where FN appears
N = size(tp,2);
figure;scatter(1:N,(fn),'filled','k')
title('False-Negative from each picture','FontSize',12)
xlabel('Image number','FontSize',12)
xticks(round(linspace(1,N,round(N/2))))
ylabel('Number of False-Negative','FontSize',12)
yticks(round(linspace(1,max(fn),max(fn))))
ylim([0 max(fn)])
s = what;
saveas(gcf,fullfile(s.path, 'False-Negative'),'fig')
saveas(gcf,fullfile(s.path, 'False-Negative'),'tif')

%%%%%%%%%%%%%%%%%% figure of where FP appears
figure;scatter(1:N,(fp),'filled','k')
title('False-Positive from each picture vs ground truth','FontSize',12)
xlabel('Image number','FontSize',12)
xticks(round(linspace(1,N,round(N/2))))
ylabel('Number of False-Positive','FontSize',12)
ylim([0 max(fn)])
yticks(round(linspace(1,max(fn),max(fn))))

saveas(gcf,fullfile(s.path, 'False-Positive'),'fig')
saveas(gcf,fullfile(s.path, 'False-Positive'),'tif')
%% %%%%%%%%%%%%%%%%%%%%% Filter the first 3 and last 3 points
% Remove data from  the first 3 and last 3 images of the whole stack
TP_filt =tp;
FN_filt = fn;
FP_filt = fp;
TP_filt(:,1:3) =[];
TP_filt(:,end-2:end) =[];
FN_filt(:,1:3) =[];
FN_filt(:,end-2:end) =[];
FP_filt(:,1:3) =[];
FP_filt(:,end-2:end) =[];
TP_filt=sum(sum(TP_filt));
FN_filt=sum(sum(FN_filt));
FP_filt=sum(sum(FP_filt));
sensitivity_filt = TP_filt/(TP_filt+FN_filt);
precision_filt = TP_filt/(TP_filt+FP_filt);
F1_filt = 2/((sensitivity_filt)^-1+(precision_filt)^-1);

%% %%%%%%%%%%%%%%%%%%%%% Save results
results.tp=tp;
results.fn=fn;
results.fp=fp;
results.TP=TP;
results.FN=FN;
results.FP=FP;
results.sensitivity=sensitivity;
results.precision=precision;
results.TP_filt=TP_filt;
results.FN_filt=FN_filt;
results.FP_filt=FP_filt;
results.sensitivity_filt=sensitivity_filt;
results.precision_filt=precision_filt;
results.F1=F1;
results.F1_filt=F1_filt;
save('results.mat', 'results');
end

function [TP,FN,idxMarkOut,valCenMark] = sensitivityFnc(I,I_ori,idxCenMark,blobMark,foldDataPred)
% Estimate then sensitivity of the centroids position on top of the images
% Sintax:
%     saveVariablesFnc(I_pred,idxCenMark,I_ori,blobMark,foldDataPred,I_mark,idxCenPred,blobPred,foldDataFP)
% Inputs:
%     I,                    Input image (predicted image stack)
%     I_ori,                Original grayscale images
%     idxCenMark,           Index of the manuel detected centrois 
%     blobMark,             Table with information of manuel detected objects
%     foldDataPred,         Folder for predicted results

% Outputs:
%     TP,                   True positive count for each section
%     FN,                   False Negative count for each section (missed
%                           data)
%     idxMarkOut,           Index of the centroids position from the Table
%     valCenMark,           Detected which centroids are corrected detected
%% Detect the position of 1
TP = zeros(1,size(I,3)); % Generate variables
FN = zeros(1,size(I,3)); % Generate variables
idxMarkOut = [];         % Generate variables
valCenMark =logical([]); % Generate variables
warning('off','all')
for i = 1:size(I,3)
    of = figure('visible', 'off'); %Overlap figure
    idxMark=find(idxCenMark(:,i)==1);
    totCen=size(idxMark,1); %total number of centroids
    overlay = imoverlay(I_ori(:,:,i),I(:,:,i), [.1 .8 .1]);
    imshow(overlay)
    axis image;
    valCen = zeros(totCen,1);
    hold on
    for j =1:totCen
        % Detect the centroid coordinate from the index
        posCen=blobMark.CentroidRound(idxMark(j),:) ;
        
        %     Detect if the value in this position is 0 or 1
        valCen(j) =logical(I(posCen(2),posCen(1),posCen(3)));
        
        hold on
        if valCen(j)==1
            scatter(posCen(1),posCen(2),'b','LineWidth',1);
        else
            scatter(posCen(1),posCen(2),'r','LineWidth',1);
        end
    end
    
    TP(i)=sum(valCen); % valCen are the centroids that are correct
    FN(i) = (totCen-sum(valCen)); % we are checking all centroid and if it equals 0 => bad
    idxMarkOut =[idxMarkOut;idxMark];
    valCenMark = [valCenMark;valCen];
    rmpath('folderthatisnotonpath')
    disp([num2str(i),' image is saved'])
    title(['Image number ',num2str(i),'.  TP: ',num2str(TP(i)),' out of ',num2str(totCen)],'FontSize', 16)
    saveas(of,[foldDataPred,'/',num2str(i),'.png']); % Save the
end
end
function [TP_pred,FP,idxMarkOut,valCenMark] = precisionFnc(I,I_ori,idxCenPred,blobPred,foldDataFP)
% Estimate then precision of the centroids position on top of the images
% Sintax:
%     [TP_pred,FP,idxMarkOut,valCenMark] = precisionFnc(I,I_ori,idxCenPred,blobPred,foldDataFP)
% Inputs:
%     I,                    Input image (predicted image stack)
%     I_ori,                Original grayscale images
%     idxCenPred,           Index of the predicted detected centrois 
%     blobPred,             Table with information of predicted detected objects
%     foldDataFP,           Folder for False positive results

% Outputs:
%     TP_pred,              True positive count for each section of the
%                           predicted images
%     FP,                   False-Positive count for each section (missed
%                           data)
%     idxMarkOut,           Index of the centroids position from the Table
%     valCenMark,           Detected which centroids are corrected detected

%% Detect the position of 1
TP_pred = zeros(1,size(I,3));
FP = zeros(1,size(I,3));
idxMarkOut = [];
valCenMark =logical([]);

warning('off','all')
for i = 1:size(I,3)
    of = figure('visible', 'off'); %Overlap figure
    idxMark=find(idxCenPred(:,i)==1);
    totCen=size(idxMark,1); %total number of centroids on this images
    overlay = imoverlay(I_ori(:,:,i),I(:,:,i), [.1 .8 .1]);
    imshow(overlay)
    axis image;
    valCen = zeros(totCen,1);
    hold on
    for j =1:totCen
        % Detect the centroid coordinate from the index
        posCen=blobPred.CentroidRound(idxMark(j),:) ;
        
        %     Detect if the value in this position is 0 or 1
        valCen(j) =logical(I(posCen(2),posCen(1),posCen(3)));
        
        hold on
        if valCen(j)==1
            scatter(posCen(1),posCen(2),'b','LineWidth',1);
        else
            scatter(posCen(1),posCen(2),'r','LineWidth',1);
        end
    end
    
    TP_pred(i)=sum(valCen); % valCen are the centroids that are correct
    FP(i) = (totCen-sum(valCen)); % we are checking all centroid and if it equals 0 => bad
    idxMarkOut =[idxMarkOut;idxMark];
    valCenMark = [valCenMark;valCen];
    rmpath('folderthatisnotonpath')
    disp([num2str(i),' image is saved'])
    title(['Image number ',num2str(i),'.  TP: ',num2str(TP_pred(i)),' out of ',num2str(totCen)],'FontSize', 16)
    saveas(of,[foldDataFP,'/',num2str(i),'.png']); % Save the
end
valCenMark = logical(valCenMark);
end



folderMain ='ValidationData_combined';
folder1 ='GrayImg';
folder2 ='ManuelMark';
folder3 ='Predicted';
resultFolder ='Results';
FN_Folder ='Results_FN';
FP_Folder ='Results_FP';

fileName = 'ori';
fileName2 = 'pred';
ImgType = 'png';

%% %%%%%%%%%%%%% Read number of files in the folder
% Number of images in the folder for analysis
N = size(dir(fullfile([folderMain,'/',folder1,'/','/*.',ImgType])),1);  
for i =1:N
    %%
    [resultPath]=folderGenerateFnc(folderMain,resultFolder);
    [FN_Path]=folderGenerateFnc(folderMain,FN_Folder);
    [FP_Path]=folderGenerateFnc(folderMain,FP_Folder);
    [I_gray] = readGrayImg(folderMain,folder1,fileName,ImgType,i);
    [I_mark] = readImg(folderMain,folder2,fileName,ImgType,i);
    [I_pred] = readImg(folderMain,folder3,fileName2,ImgType,i);
    % Save quantitative measurement in the results folder
    Results =resultEstimationFnc(I_mark,I_pred,resultPath,i);
    % Save figures
    saveFigureFnc(I_mark,I_pred,resultPath,FP_Path,FN_Path,i)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [folderPath]=folderGenerateFnc(folderMain,folderName)
% Create the save folder of images
% Sintax:
%     [folderSave]=folderGenerateFnc(folder1,folder2)
% Inputs:
%     folder1,     First folder name
%     name,        Name of image file

% Outputs:
%     folderSave,         Folder name for saved images

%     srcFiles = dir(fullfile(folderMain));
%     path=srcFiles.folder;
%     path=[path '\'];
    
s = what(folderMain);
savePath=s.path;
folderPath = [savePath,'\',folderName];
if ~exist(folderPath, 'dir')
    mkdir(folderPath) % Creating a folder with the brain ID
end


    
end

function [I] = readImg(folderMain,folder,fileName,ImgType,i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
srcFiles = dir([folderMain,'/',folder,'/*.',ImgType]);  % the folder in which ur images exists
% read the first image
path=srcFiles.folder;
path=[path '\'];
% Read image
I = imread([path,fileName,' (',num2str(i),').',ImgType]);
%% Convert image to grayimage and load images in folder
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[~, ~, numberOfColorChannels] = size(I);
if numberOfColorChannels > 1
    % Convert RGB to gray
    I = imbinarize(rgb2gray((I)));
else
    I=logical(I);
end

%% Remove rows and columns of uneven dimension
if mod(size(I,1),2) ==1
    I(end,:,:)=[];
end
if mod(size(I,2),2) ==1
    I(:,end,:) =[];
end
disp('Binary image loaded')
% I = imfill(I,'holes');
end


function [I_gray] = readGrayImg(folderMain,folder1,fileName,ImgType,i)
srcFiles = dir([folderMain,'/',folder1,'/*.',ImgType]);  % the folder in which ur images exists
% read the first image
path=srcFiles.folder;
path=[path '\'];
% Read image
I = imread([path,fileName,' (',num2str(i),').',ImgType]);
%% Convert image to grayimage and load images in folder
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[~, ~, numberOfColorChannels] = size(I);
if numberOfColorChannels > 1
    % Convert RGB to gray
    I_gray = (rgb2gray((I)));
else
    I_gray=(I);
end
disp('Gray image loaded ')

end


function Results =resultEstimationFnc(I_mark,I_pred,folderPath,i)
Results.TP =sum(sum(I_mark.*I_pred));
    Results.TN =sum(sum((1-I_mark).*(1-I_pred)));
    Results.FP =sum(sum((1-I_mark).*I_pred));
    Results.FN =sum(sum(I_mark.*(1-I_pred)));
    Results.Sensitivity = Results.TP/(Results.TP+Results.FN);
    Results.Precision = Results.TP/(Results.TP+Results.FP);
    Results.F1 = 2/((Results.Sensitivity)^-1+(Results.Precision)^-1);
    Results.Accuracy = (Results.TP+Results.TN)/(Results.TP+Results.TN+Results.FP+Results.FN);
    save(fullfile(folderPath, ['Results',num2str(i),'.mat']), 'Results')
    
end



function saveFigureFnc(I_mark,I_pred,resultPath,FP_Path,FN_Path,i)
 %% Save figures
    % Enlarge figure to full screen.
    f = figure('visible', 'off');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(1,4,1)
    imshow(I_mark)
    title('Mark')
    subplot(1,4,2)
    imshow((I_pred))
    title('Pred')
    subplot(1,4,3)
    imshow((1-I_mark).*I_pred)
    title('FP')
    subplot(1,4,4)
    imshow(I_mark.*(1-I_pred))
    title('FN')
    saveas(gcf,fullfile(resultPath, ['Brain',num2str(i)]),'fig')
    saveas(gcf,fullfile(resultPath, ['Brain',num2str(i)]),'png')
    close(f)

  
    
     f = figure('visible', 'off');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow((1-I_mark).*I_pred)
    saveas(gcf,fullfile(FP_Path, ['FP',num2str(i)]),'fig')
    saveas(gcf,fullfile(FP_Path, ['FP',num2str(i)]),'png')
     close(f)
     
         f = figure('visible', 'off');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow(I_mark.*(1-I_pred))
    saveas(gcf,fullfile(FN_Path, ['FN',num2str(i)]),'fig')
    saveas(gcf,fullfile(FN_Path, ['FN',num2str(i)]),'png')
    close(f)
end



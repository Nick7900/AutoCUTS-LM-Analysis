% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read files %%%%%%%%%%%%%%%%%%%%%%%%%%%
MANUEL = 1; % Manuel crop image
name = 'Example';
folder1 ='AutuCUTS_Pipeline';
folder2 =[name,'_2_align'];
fileName = 'Img_align_0';
ImgType = 'png';
ImgAmount = 3;
buffer = 0;

% Define the cropping size of the images
% Increase or decrease with 256 pixels for the width and height
factor = -1; % Define the factor to increase or decrease of ROI
%%%%%%%% Variables
height =(4096)+256*factor;
width=(4096)+256*factor;


%% Start program %%%%%%%%%%%%%%%
% Determine where to crop images
imgNum = 2; %read image number out of stack
[I,N] = readImgSingleFnc(folder1,folder2,fileName,ImgType,imgNum);
% [I,N] = readMultipleImagesFnc(folder1,folder2,fileName,ImgType,ImgAmount,buffer);
%%%%%%%%%% Choose top left corner from cropping %%%%%%%%%%%%%
[xmin,ymin] =visCropFncFnc(I,height,width,MANUEL); % visualize crop image area
pause(0.1)
% [xmin,ymin] =visCropManyFnc(I,height,width,MANUEL);
cRect=[xmin(1) ymin(1) width height]; % ROI of the crop image
%% Crop all images and save
%%%%%%%%%% Save Images
[fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1,cRect);

% Save cropped images
for i = 1:N
    x =readImgSingleLoopFnc(folder1,folder2,fileName,ImgType,cRect,N,i+buffer,buffer);
    imwrite(x, [folderNameSave,'\',fileNameSave,' (',num2str(i+buffer),')','.png'],'png');
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,N] = readImgSingleFnc(folder1,folder2,fileName,ImgType,imgNum)
% Read one image to determine cropping area
% Sintax:
%     [I,N] = readImgSingleFnc(folder1,folder2,fileName,ImgType,imgNum)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     fileName,    Name of image that we read from
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     imgNum,      Read image number X of the image stack

% Outputs:
%     I,            Output image name
%     N,            Number of images in the folder

srcFiles = dir([fullfile(folder1,'/',folder2),'/*.',ImgType]);
path=srcFiles.folder;
path=[path '\'];
N =size(srcFiles,1);
I_mark = imread([path,fileName,num2str(imgNum),'.',ImgType]);
%  Convert image to grayimage and load images in folder
[~, ~, numberOfColorChannels] = size(I_mark);
if numberOfColorChannels > 1
    disp(['Visualize image number ',num2str(imgNum),' out of ',num2str(N)])
    I = rgb2gray(imread([path,fileName,num2str(imgNum),'.',ImgType]));
else
    disp(['Visualize image number ',num2str(imgNum),' out of ',num2str(N)])
    I = (imread([path,fileName,num2str(imgNum),'.',ImgType]));
    
end

end


function [I,N] = readMultipleImagesFnc(folder1,folder2,fileName,ImgType,ImgAmount,buffer)
% Read one image to determine cropping area
% Sintax:
%     [I,N] = readMultipleImagesFnc(folder1,folder2,fileName,ImgType,ImgAmount,buffer)
%
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     fileName,    Name of image that we read from
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     ImgAmount,   How many images of the total stack of images are used
%     buffer,      Buffer to change the starting number of images

% Outputs:
%     I,            Output image name
%     N,            Number of images in the folder

% Get path of files
srcFiles = dir([fullfile(folder1,'/',folder2),'/*.',ImgType]);
path=srcFiles.folder;
path=[path '\'];
N =size(srcFiles,1);
I_mark = imread([path,srcFiles(1).name]);

%  Convert image to grayimage and load images in folder
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[~, ~, numberOfColorChannels] = size(I_mark);
I = uint8(zeros(size(I_mark,1),size(I_mark,2),ImgAmount));
if numberOfColorChannels > 1
    ImgDivide =N/(N*ImgAmount);
    if ImgAmount==1
        disp(['Load image number ',num2str(round(N/2)),' out of ',num2str(N)]);
        I = rgb2gray(imread([path,srcFiles(round(N/2)).name]));
    else
        for i = 1:ImgAmount
            if i ==ImgAmount
                round((ImgDivide*i)*N);
                disp(['Load image number ',num2str(size(srcFiles,1)),' out of ',num2str(N)])
                I(:,:,ImgAmount) = rgb2gray(imread([path,srcFiles(end).name]));
            else
                imgNum = round((ImgDivide*i)*N);
                disp(['Load image number ',num2str(imgNum),' out of ',num2str(N)]);
                I(:,:,i) = rgb2gray(imread([path,srcFiles(imgNum).name]));
            end
        end
    end
else
    if ImgAmount==1
        disp(['Load image number ',num2str(round(N/2)),' out of ',num2str(N)]);
        I(:,:,ImgAmount) = imread([path,fileName,' (',num2str(N/2),').png']);
    else
        ImgDivide =N/(N*ImgAmount);
        for i = 1:ImgAmount
            if i ==ImgAmount
                disp(['Load image number ',num2str(size(srcFiles,1)+buffer),' out of ',num2str(N+buffer)])
                I(:,:,ImgAmount) = imread([path,fileName,num2str(N+buffer),'.png']);
            else
                imgNum = round((ImgDivide*i)*N)+buffer;
                disp(['Load image number ',num2str(imgNum),' out of ',num2str(N+buffer)]);
                I(:,:,i) = imread([path,fileName,num2str(imgNum),'.png']);
            end
        end
    end
end
end


function [fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1,cRect)
% Create the save folder of images
% Sintax:
%     [fileNameSave,folderNameSave]=folderGenerateFnc(name,folder1)
% Inputs:
%     name,        name of subject
%     folder1,     First folder name
%     cRect,       Size and position of the crop rectangle in spatial coordinates

% Outputs:
%     fileNameSave,         Output image name
%     folderNameSave,       Directory for saved folder

s = what(folder1);
savePath=s.path;
folderSave =[name,'_3_crop_DeepLearning_x',num2str(cRect(1)),'_y',num2str(cRect(2)),];
folder = [savePath,'\',folderSave];
fileNameSave ='ori';
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end
%%%%%%%%%%%% Save coordinate crop positions
save(fullfile(savePath, 'cRectDeep.mat'), 'cRect')

end


function [x] = readImgSingleLoopFnc(folder1,folder2,fileName,ImgType,cRect,N,i,buffer)
% Read one image of the time from folder and crop it
% Note - The computer was out of memory but reading and cropping the 
% whole stack, which is the reason for this approach. 
% Sintax:
%     [x] = readImgSingleLoopFnc(folder1,folder2,fileName,ImgType,cRect,N,i,buffer)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     fileName,    Name of image that we read from
%     ImgType,     Image type (eg. TIF, png, JPEG)
%     cRect,       Size and position of the crop rectangle in spatial coordinates
%     N,           Number of images in the folder
%     i,            Image number from loop
%     buffer,      Buffer to change the starting number of images

% Outputs:
%     x,                    Output image name

srcFiles = dir([folder1,'/',folder2,'/*.',ImgType]);  % the folder in which ur images exists
path=srcFiles.folder;
path=[path '\'];
disp(['Load image number ',num2str(i),' out of ',num2str(N+buffer)])
I = (imread([path,fileName,num2str(i),'.',ImgType]));

%% Crop Image
x = cropImgLoopFnc(I,cRect);
end

function [x,y] =visCropFncFnc(I,height,width,MANUEL)
% Define cropping area
% Sintax:
%     [x,y] =visCropFncFnc(I,height,width,MANUEL)
% Inputs:
%     I,            Image used to choose cropping area
%     height,       Define height of cropping area
%     width,    	Define width of cropping area
%     MANUEL,       Define if the user manuel is going to choose cropping
%                   area

% Outputs:
%     x,            Size and position of the crop rectangle in x-axis
%     y,            Size and position of the crop rectangle in y-axis

if MANUEL ==1
    imshow(I)
    title('Press at top left corner to visualize cropping ')
    lineWidth = 2;
    [xi, yi] = ginput(1);
    % Visualize crop output
    % p1 = [xi; yi];
    % p2 = [xi; yi+height];
    % p3 = [xi+width; yi];
    % p4 = [xi+width; yi+height];
    x = [xi xi xi+width xi+width];
    y = [yi yi+height yi yi+height];
    hold on

    pause(0.5)
    plot([xi xi+width],[yi yi],'--y','LineWidth',lineWidth) % top line
    plot([xi xi],[yi yi+height],'--y','LineWidth',lineWidth) % left line
    plot([xi xi+width],[yi+height yi+height],'--y','LineWidth',lineWidth) % bottom line
    plot([xi+width xi+width],[yi yi+height],'--y','LineWidth',lineWidth) % right line
    scatter(x,y,'r','filled','LineWidth',lineWidth)
else
    xi = 1;
    yi =1;
    x = [xi xi xi+width xi+width];
    y = [yi yi+height yi yi+height];
    
end
end

function [xi_out,yi_out] =visCropManyFnc(I,height,width,MANUEL)
% Define cropping area
% Sintax:
%     [xi_out,yi_out] =visCropManyFnc(I,height,width,MANUEL)
% Inputs:
%     I,            Image used to choose cropping area
%     height,       Define height of cropping area
%     width,    	Define width of cropping area
%     MANUEL,       Define if the user manuel is going to choose cropping
%                   area

% Outputs:
%     xi_out,            Size and position of the crop rectangle in x-axis
%     yi_out,            Size and position of the crop rectangle in y-axis
if MANUEL ==1
    xi = zeros(size(I,3),1);
    yi = zeros(size(I,3),1);
    for i = 1:size(I,3)
        imshow(I(:,:,i))
        [xi(i), yi(i)] = ginput(1);
        % Visualize crop output
        % p1 = [xi; yi];
        % p2 = [xi; yi+height];
        % p3 = [xi+width; yi];
        % p4 = [xi+width; yi+height];
        hold on
        x = [xi(i) xi(i) xi(i)+width xi(i)+width];
        y = [yi(i) yi(i)+height yi(i) yi(i)+height];
        scatter(x,y,'r','filled')
        pause(2)
    end
    xi_out = round(mean(xi));
    yi_out = round(mean(yi));
    
    if size(I,3)<=2
        figure(2)
        subplot(1,2,i)
        imshow(I(:,:,i))
        hold on
        x = [xi_out xi_out xi_out+width xi_out+width];
        y = [yi_out yi_out+height yi_out yi_out+height];
        scatter(x,y,'r','filled')
        
    else
        for i = 1:size(I,3)
            figure(2)
            subplot(2,2,i)
            imshow(I(:,:,i))
            hold on
            x = [xi_out xi_out xi_out+width xi_out+width];
            y = [yi_out yi_out+height yi_out yi_out+height];
            scatter(x,y,'r','filled')
        end
    end
else
    x = 1;
    y =1;
    xi_out = [x x x+width x+width];
    yi_out = [y y+height y y+height];
    
end
end

function x = cropImgLoopFnc(I,cRect)
% Crop image
% Sintax:
%     x = cropImgLoopFnc(I,cRect)
% Inputs:
%     I,           Image input that will be cropped
%     cRect,       Size and position of the crop rectangle in spatial coordinates

% Outputs:
%     x,           Cropped output image

cRect = round(cRect);
x = imcrop(I,cRect);
%% Remove rows and columns
if mod(size(x,1),2) ==1
    x(end,:,:)=[];
end
if mod(size(x,2),2) ==1
    x(:,end,:) =[];
end

end

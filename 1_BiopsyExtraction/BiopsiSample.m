% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% clc; clear;
LINEMANUEL = false;
name='ImgExample';
folder1 = 'Analyse';
ImgType ='png';
N_points = 1; % Number of sample points for each quartile

% Reading the image and directory information
[I,name,srcFiles,path]=imgReadInfoFnc(name,folder1,ImgType);

% Generate folder to save files
[folderSave]=folderGenerateFnc(folder1,name);

% Define biopsy width on tissue
[radii] = getRadiiValue(I,LINEMANUEL);

% Define region of interest for the analysis
[I_gray,I_sample] =cropAreaFnc(I);
%%%%%%%%%%%%%%%%%%%%% Run analysis
[I_EroBW,C,eroImg,outPoint1,outPoint2] = sampleFncUpdate(I_sample, N_points, radii, name,srcFiles); % Output

saveData(I,I_gray,I_sample,eroImg,C,outPoint1,outPoint2,I_EroBW,radii,path,name)

%% %%%%%%%%%%%%%%%%%%%%% Functions
function [folderSave]=folderGenerateFnc(folder1,name)
% Create the save folder of images
% Sintax:
%     [folderSave]=folderGenerateFnc(folder1,folder2)
% Inputs:
%     folder1,     First folder name
%     name,        Name of image file

% Outputs:
%     folderSave,         Folder name for saved images

s = what(folder1);
savePath=s.path;
folderSave =folder1;
folderName = [savePath,'\',folderSave];
if ~exist(folderName, 'dir')
    mkdir(folder1,name) % Creating a folder with the brain ID
end
end


function [I,name,srcFiles,path]=imgReadInfoFnc(name,folder1,ImgType)
% Read image and directory information
% Sintax:
%     [folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType)
% Inputs:
%     name,        Name of image file
%     folder1,     First folder name
%     ImgType,     Image type (eg. TIF, png, JPEG)

% Outputs:
%     I,           Image output
%     name,        Name of image file
%     srcFiles,    Structure with file informations
%     path,        Directory of the image position
try
    srcFiles = dir(fullfile(folder1));
    path=srcFiles.folder;
    path=[path '\'];
    I =imread(strcat(path, name,'.png'));
catch
    [file,path] = uigetfile(['*.',ImgType]);
    name =erase(file,['.',ImgType]);
    if isequal(file,0)
        disp('User selected Cancel');
    else
        disp(['User selected ', fullfile(path,file)]);
    end
    % Read images
    I =imread(strcat(path, name,'.png'));
end

end

function radii = getRadiiValue(I,LINEMANUEL)
% Define length of biopsi size to define sample area
% Sintax:
%     radii = getRadiiValue(I,LINEMANUEL)
% Inputs:
%     I,             Image input
%     LINEMANUEL,    Using default value or manually choose a new value

% Outputs:
%     radii,         Radius of the biopsy

if LINEMANUEL
    captionFontSize = 16;
    figure, imshow(I);title('Press enter to save distance value','FontSize',captionFontSize)
    h = imdistline(gca);
    api = iptgetapi(h);
    % wait for key pressed
    pause('on');
    pause;
    % get line
    distance = api.getDistance();
    radii = distance/2;
    close all
else
    % By default a width of 5 mm got a distance of 60 pixels
    diameter = 60;
    radii = diameter/2;
end
end


function [I_gray,I_sample] =cropAreaFnc(I)
% Define region of interest to sample the biopsies
% Sintax:
%     [I_gray,I_thresh] =cropAreaFnc(I)
% Inputs:
%     I,             Image input

% Outputs:
%     I_gray,         Grayscale image of the input image
%     I_sample,       Sample area of the input image
%% Removing the distance line
captionFontSize = 16;
I_gray = rgb2gray(I);
imshow(I_gray);
title('Specify region of interest','FontSize',captionFontSize)
h = drawpolygon('FaceAlpha',0); %draw something
% true/1/white inside the drawn area and false/0/black outside the drawn area.
M = ~h.createMask(); % ~ inverts that mask so it's black inside the drawn area and white outside.
I_crop= I_gray;
I_crop(M) = 0;
imshow(I_crop);
title('Region of interest','FontSize',captionFontSize)

%% Select boundaries to flat areas
p = 2; % Chosen points
title(['Mark ',num2str(p),' darkest pixels near the flat surface'],'FontSize',16 )
[py,px] = ginput(p);
py = round(py);
px = round(px);
index =zeros(1,p);
for i =1:1:p
    index(i) =double(I_crop(px(i,:),py(i,:)));
end
pause('on');
pause('off');
close all
avgIndex =mean(index);

I_sample = I_crop;
% Set values below boundary to 0
I_sample(I_sample<=avgIndex)=0;
end




function [I_EroBW,C,eroImg,outPoint1,outPoint2] = sampleFncUpdate(I_sample, N_points, radii, name,srcFiles)
% Random choose sample points
% Sintax:
%     [I_EroBW,C,eroImg,outPoint1,outPoint2] = sampleFncUpdate(I_sample, N_points, radii, name,srcFiles)
% Inputs:
%     I_sample,       Sample area of the input image
%     N_points,       Number of sample points for each quartile
%     radii,          Radius of the biopsy
%     name,           Name of image file
%     srcFiles,       Structure with file informations

% Outputs:
%     I_EroBW,        Binary sample area
%     C,              Combined image of the sample area+grayscale image
%     eroImg,         Erosion image in grayscale
%     outPoint1,      Points for the biopsies in the first quartile
%     outPoint2,      Points for the biopsies in the second quartile

captionFontSize = 16; % Size of text of plots
%% Erosin
SE = offsetstrel('ball',round(radii),round(radii)); % define the filter size (ball)
eroImg = imerode(I_sample,SE);
% Making a threshold and make it black
% thres =max(max(IM1))*0.66;
% IM1(IM1<thres)=0;
% Finding the threshold for the binary image
% Binary images
level = graythresh(eroImg);
I_EroBW = imbinarize(eroImg,level);
figure(1);
imshow(I_EroBW)

fname = srcFiles.folder;
directory_out = [fname,'\',name];
str = '_BWImg';
OutputName = strcat(name,str);
saveas(gca, fullfile(directory_out,OutputName),'fig');
saveas(gca, fullfile(directory_out,OutputName),'png');
% Fuse images
Combined = imfuse(I_sample,eroImg,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
S1 = [max(size(Combined,2))*0.75 max(size(Combined,1))*0.78];
S2 = [max(size(Combined,2))*0.75 max(size(Combined,1))*0.88];

position = [S1;S2];
box_color = {'yellow','red'};
text_str = cell(2,1);
text_str{2} = 'Filtered Area';
text_str{1} = 'Sample Area';

C = insertText(Combined,position,text_str,'FontSize',25,'BoxColor',...
    box_color,'BoxOpacity',0.4,'TextColor','white');
figure(2);
imshow(C);
title('Fused image', 'FontSize', captionFontSize)
hold off;

str = '_Fused';
OutputName = strcat(name,str);
% Saving the data
saveas(gca, fullfile(directory_out,OutputName),'fig');
saveas(gca, fullfile(directory_out,OutputName),'png');
hold off
%% Analysing image
[sX, sY]=size(I_EroBW);  % Getting the size of the image
% sumRow =[]
sumRow=zeros(1,sX); % Empty matrix that will get the sum value for each row
for k=1:1:sX
    sumRow(1,k)=sum(I_EroBW(k,:)); % row number and sum all columns
end

cumRow=cumsum(sumRow);
% Sum all white pixels
nWhite=sum(sumRow);    % total amount white pixels

% Find first white pixel
ny0 =(find(sumRow > 0, 1));

% Remember Y is row hight
[~,ny25]=find(cumRow==max(cumRow(cumRow<=nWhite*.25))); % Marking the row that contains 25 %
[~,ny50]=find(cumRow==max(cumRow(cumRow<=nWhite*.50))); % Marking the row that contains 50 %
[~,ny75]=find(cumRow==max(cumRow(cumRow<=nWhite*.75))); % Marking the row that contains 75 %
[ny100]=find(cumRow == max(cumRow(:)),1,'first'); % Marking the row that contains 75 %

%% Defining stripes
nWhite=sum(sumRow);    % total amount white pixels
% Remember Y is row hight
[~,nx25]=find(cumRow==max(cumRow(cumRow<=nWhite*.25))); % Marking the row that contains 25 %
[~,nx50]=find(cumRow==max(cumRow(cumRow<=nWhite*.50))); % Marking the row that contains 50 %
[~,nx75]=find(cumRow==max(cumRow(cumRow<=nWhite*.75))); % Marking the row that contains 75 %

% Dividing the image in 4 quartiles
x25 = [nx25 nx25];
x50 = [nx50 nx50];
x75 = [nx75 nx75];
y = [0 max(cumRow)]; %higest value of a row

figure(3);
caption = sprintf('White pixel count');
title(caption, 'FontSize', captionFontSize);
yyaxis left
plot(1:1:sX,sumRow)
ylabel('Sum white pixels', 'FontSize',10, 'FontWeight', 'Bold') % left y-axis
xlabel('Row number', 'FontSize',13, 'FontWeight', 'Bold')
yyaxis right
plot(1:1:sX,cumRow)
ylabel('Cumulative pixels', 'FontSize',10, 'FontWeight', 'Bold') % right y-axis

hold on % new figure
xt = [nx25 nx50 nx75];
yt = [max(cumRow) max(cumRow) max(cumRow)];
str = {'25%','50%','75%'};
text(xt,yt,str)

line(x25,y,'Color','g','LineStyle','--')
line(x50,y,'Color','g','LineStyle','--')
line(x75,y,'Color','g','LineStyle','--')
hold off;

%% Checking acquired points correspond to all white points
[nyWhite,nxWhite]=find(I_EroBW>0);
%% Finding all white points at different quarters
[row25]=find(nyWhite<ny25);         % points in 1st quarter area(rows)
[row50_up]=find(nyWhite<ny50);          % points in 2nd quarter area
[row50_down]=find(nyWhite>ny25);
row50=intersect(row50_up,row50_down);
[row75_up]=find(nyWhite<ny75);          % points in 3rd quarter area
[row75_down]=find(nyWhite>ny50);
row75=intersect(row75_up,row75_down);
[row100]=find(nyWhite>ny75);            % points in 4th quarter area

P25=[nxWhite(row25) nyWhite(row25)];        % 1st quarter, x and y
P50=[nxWhite(row50) nyWhite(row50)];        % 2nd quarter, x and y
P75=[nxWhite(row75) nyWhite(row75)];        % 3rd quarter, x and y
P100=[nxWhite(row100) nyWhite(row100)];     % 4th quarter, x and y

nP25=(1:1:length(P25));             % scatter points across 1st quarter area sector
nP50=(1:1:length(P50));             % scatter points across 2nd quarter area sector
nP75=(1:1:length(P75));             % scatter points across 3rd quarter area sector
nP100=(1:1:length(P100));           % scatter points across 4th quarter area sector

str = '_Histogram';
OutputName = strcat(name,str);
% Saving the data
saveas(gca, fullfile(directory_out,OutputName),'fig');
saveas(gca, fullfile(directory_out,OutputName),'png');
%% Sampling of the quartiles
sampleR = randi([0 1]);
% Applying the quartiles to the image
%Coordinate of the image
[~, top_row]= find(sum(I_EroBW,3)', 1);
[~, bottom_row]= find(sum(I_EroBW,3)', 1, 'last');

% getting the dimension of the picture
[~,nI] = size(I_EroBW);

% Plotting the Patches
p25 = [0 nx25; nI nx25; nI top_row; 0 top_row];
p50 = [0 nx50; nI nx50; nI nx25; 0 nx25];
p75 = [0 nx75; nI nx75; nI nx50; 0 nx50];
p100 = [0 bottom_row; nI bottom_row; nI nx75; 0 nx75];

if sampleR ==1
    figure(4);
    imshow(I_EroBW);
    title('Sample area', 'FontSize', captionFontSize);hold on;
    s=1;
    str_color=['b' 'r' 'y' 'r' 'm' 'c' 'w' 'k'];
    str_marker=['o' 'd' 's' '+' '*' '.' '^' 'v' 'p' 'h'];
    Psct25=[];  % counter for scattered points across 1st quarter area
    Psct75=[];  % counter for scattered points across 3rd quarter area
    hold on
    % patch(vertices(:,1), vertices(:,2), 'r', 'EdgeAlpha', 0.1);
    h25 = patch(p25(:,1), p25(:,2), 'r');
    h50 = patch(p50(:,1), p50(:,2), 'white');
    h75 = patch(p75(:,1), p75(:,2), 'r');
    h100 = patch(p100(:,1), p100(:,2), 'white');
    
    hold on;
    alpha(h25,0.2)
    alpha(h50,0.5)
    alpha(h75,0.2)
    alpha(h100,0.5)
    % saveas(gcf,'quartile.fig')
    hold on;
    
    plot(ones(1,sY)*ny0,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny25,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny50,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny75,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny100,'r', 'LineWidth',2);
    plot(P25(:,1),P25(:,2),'Color',[1 .5 .5],'LineStyle','none','Marker','.')   % check points 1st quarter area
    plot(P75(:,1),P75(:,2),'Color',[1 .5 .5],'LineStyle','none','Marker','.')   % check points 3rd quarter area
    
    while s<=N_points
        n25=randi(numel(nP25),1,1);
        n75=randi(numel(nP25),1,1);
        hold all;
        plot(P25(n25,1),P25(n25,2),'LineStyle','none','Color',str_color(s),'Marker',str_marker(s), 'LineWidth',2);   % (d)iamond o:circle (s)quare + * x . ^ v < > (p)entagram (h(exagram
        plot(P75(n75,1),P75(n75,2),'LineStyle','none','Color',str_color(s),'Marker',str_marker(s), 'LineWidth',2);
        Psct25=[Psct25;P25(n25,:)];   % recording scattered points
        Psct75=[Psct75;P75(n75,:)];   % recording scattered points
        P25(n25,:)=[];nP25(n25)=[];
        P75(n75,:)=[];nP75(n75)=[];
        s=s+1;
    end
    hold off
    outPoint1 =Psct25;
    outPoint2 =Psct75;
    str = '_SamplePoints';
    f1 = 'Analyse/';
    f2 = '/';
    newName = strcat(name,f2);
    filename = strcat(f1,newName);
    OutputName = strcat(name,str);
    % Saving the data
    saveas(gca, fullfile(filename,OutputName),'fig');
else
    figure(4);
    imshow(I_EroBW);
    title('Sample area', 'FontSize', captionFontSize);hold on;
    s=1;
    str_color=['b' 'r' 'y' 'r' 'm' 'c' 'w' 'k'];
    str_marker=['o' 'd' 's' '+' '*' '.' '^' 'v' 'p' 'h'];
    Psct50=[zeros(1,2)];  % counter for scattered points across 2nd quarter area
    Psct100=[zeros(1,2)];  % counter for scattered points across 4th quarter area
    
    h25 = patch(p25(:,1), p25(:,2), 'white');
    h50 = patch(p50(:,1), p50(:,2), 'r');
    h75 = patch(p75(:,1), p75(:,2), 'white');
    h100 = patch(p100(:,1), p100(:,2), 'r');
    
    hold on;
    alpha(h25,0.5)
    alpha(h50,0.2)
    alpha(h75,0.5)
    alpha(h100,0.2)
    % saveas(gcf,'quartile.fig')
    hold on;
    
    
    plot(ones(1,sY)*ny0,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny25,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny50,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny75,'r', 'LineWidth',2);
    plot(ones(1,sY)*ny100,'r', 'LineWidth',2);
    plot(P50(:,1),P50(:,2),'Color',[1 .5 .5],'LineStyle','none','Marker','.')   % check points 1st quarter area
    plot(P100(:,1),P100(:,2),'Color',[1 .5 .5],'LineStyle','none','Marker','.')   % check points 1st quarter area
    
    while s<=N_points
        n50=randi(numel(nP50),1,1);
        n100=randi(numel(nP50),1,1);
        hold all;
        plot(P50(n50,1),P50(n50,2),'LineStyle','none','Color',str_color(s),'Marker',str_marker(s), 'LineWidth',2);   % (d)iamond o:circle (s)quare + * x . ^ v < > (p)entagram (h(exagram
        plot(P100(n100,1),P100(n100,2),'LineStyle','none','Color',str_color(s),'Marker',str_marker(s), 'LineWidth',2);
        Psct50=[Psct50;P50(n50,:)];   % recording scattered points
        Psct100=[Psct100;P100(n100,:)];   % recording scattered points
        P50(n50,:)=[];nP50(n50)=[];
        P100(n100,:)=[];nP100(n100)=[];
        s=s+1;
        
    end
    hold off;
    outPoint1 =Psct50;
    outPoint2 =Psct100;
    % saveas(gcf,'quartile.fig')
    
    str = '_SamplePoints';
    OutputName = strcat(name,str);
    % Saving the data
    saveas(gca, fullfile(directory_out,OutputName),'fig');
    saveas(gca, fullfile(directory_out,OutputName),'png');
end

end

function saveData(I,I_gray,I_sample,eroImg,C,outPoint1,outPoint2,I_EroBW,radii,path,name)
% Save all variables to the directory
% Sintax:
%     saveData(I,I_gray,I_sample,eroImg,C,outPoint1,outPoint2,I_EroBW,radii,path,name)
% Inputs:
%     I,              Input image
%     I_gray,         Grayscale image
%     I_sample,       Sample area of the input image
%     eroImg,         Erosion image in grayscale
%     C,              Combined image of the sample area+grayscale image
%     outPoint1,      Points for the biopsies in the first quartile
%     outPoint2,      Points for the biopsies in the second quartile
%     I_EroBW,        Binary sample area
%     radii,          Radius of the biopsy
%     path,        Directory of the image position
%     name,           Name of image file

% Save data for further use in the future
field1 = 'I';        value1= I;
field2 = 'I_gray';    value2 =I_gray;
field3 = 'I_sample';      value3 =I_sample;
field4 = 'ero';         value4 =eroImg;
field5 = 'C';           value5 =C;
field6 = 'point1';      value6 =outPoint1;
field7 = 'point2';      value7 =outPoint2;
field8 = 'I_EroBW';      value8 =I_EroBW;
Data = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8);

% Saving the data
str = '_Data';
OutputName = strcat(str);
save([[path,name], OutputName '.mat'], 'Data');
%%
captionFontSize = 16;
% Visualize the biopshy points
f = figure('visible', 'off');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
imshow(I)
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', 40);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', 40);
caption = sprintf('Raw image');
title(caption, 'FontSize', captionFontSize);
% Saving the data
str = '_raw';
OutputName = strcat(name,str);
saveas(gca, fullfile([path,name],OutputName),'png');
saveas(gca, fullfile([path,name],OutputName),'fig');
close(f)

% imfused image
f = figure('visible', 'off');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
imshow(C)
hold on
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
caption = sprintf('Imfused image');
title(caption, 'FontSize', 12);
% Saving the data
str = '_imfused';
OutputName = strcat(name,str);
saveas(gca, fullfile([path,name],OutputName),'png');
saveas(gca, fullfile([path,name],OutputName),'fig');
close(f)


% Original +Imfused image
f = figure('visible', 'off');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
imshow(I)
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
caption = sprintf('Raw image');
title(caption, 'FontSize', 12);
subplot(1,2,2)
imshow(C)
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
caption = sprintf('Greyscale picture(red), Binary(yellow)');
title(caption, 'FontSize', 12);
% Saving the data
str = '_comparison';
OutputName = strcat(name,str);
saveas(gca, fullfile([path,name],OutputName),'png');
saveas(gca, fullfile([path,name],OutputName),'fig');
close(f)

f = figure('visible', 'off');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
imshow(I_gray)
title('Grey image')
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
subplot(1,3,2)
imshow(I_sample)
title('Threshold image')
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
subplot(1,3,3)
imshow(eroImg)
title('Erosion')
hold on;
plot(outPoint1(:,1), outPoint1(:,2), 'b.', 'MarkerSize', radii);
plot(outPoint2(:,1), outPoint2(:,2), 'b.', 'MarkerSize', radii);
str = '_grayCombine';
OutputName = strcat(name,str);
saveas(gca, fullfile([path,name],OutputName),'png');
saveas(gca, fullfile([path,name],OutputName),'fig');
close(f)
end
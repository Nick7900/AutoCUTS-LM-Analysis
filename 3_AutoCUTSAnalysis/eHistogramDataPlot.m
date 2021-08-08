% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% Measured properties of the pyramidal neurons
voxelRes = 272; % nm
z = 800; % nm
name = 'subject1';
folder1 ='AutuCUTS_Pipeline';
folder2 ='Example_4_Prediction';
folder3 ='Example_UNetDense_stack_layer3_results';

% Define folders
[folder1,folder2,folder3]=fileInforFnc(folder1,folder2,folder3);
% Read files
[folder,blobfiltPyramid,pyraRadius2D,windowRes]=readFilesFnc(folder1,folder2,folder3,name,voxelRes);
% Plot histograms and save variables
[cellTable]=plotEstimationFnc(folder,windowRes,name,voxelRes,z,blobfiltPyramid,pyraRadius2D);

%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
function [folder1,folder2,folder3,srcFiles,path]=fileInforFnc(folder1,folder2,folder3)
%%
% Get path to read files
% Sintax:
%     [folder2,srcFiles]=fileInforFnc(folder1,folder2,ImgType)
% Inputs:
%     folder1,     First folder name
%     folder2,     Second folder name
%     folder3,     Third folder name

% Outputs:
%     srcFiles,     Structure with file informations
%     folder1,      First folder name
%     folder2,      Second folder name
%     folder3,      Third folder name
try
    srcFiles = dir([fullfile(folder1,'/',folder2,'/',folder3),'/*.']);
    path=srcFiles.folder;
    path=[path '\'];
    if isempty(srcFiles)==1
        path = uigetdir;
        srcFiles=dir([fullfile(path),'/*.']);
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
        srcFiles=dir([fullfile(path),'/*.']);
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


function [folder,blobfiltPyramid,pyraRadius2D,windowRes]=readFilesFnc(folder1,folder2,folder3,name,voxelRes)
% Read files, create save folders and determine window
% Sintax:
%     [folder,blobfiltPyramid,windowRes]=readFilesFnc(folder1,folder2,folder3,name,voxelRes)
% Inputs:
%     folder1,          First folder name
%     folder2,          Second folder name
%     folder3,          Third folder name
%     name,             name of subject/folder that will be analysed
%     voxelRes,         Voxel resolution of the analysed data

% Outputs:
%     folder,           Folder name for saved images
%     blobfiltPyramid,  Table with information of detected pyramidal cells
%     windowRes,        Defined sampling window adjusted with the voxel
%                       resolution

%%%%%%% Create save folder
s = what(folder1);
savePath=[s.path,'\',folder2];
folderSave =['Example_UNetDense_',name,'_histogram_results'];
folder = [savePath,'\',folderSave];
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end

% Load data
folderData = [savePath,'\',folder3];
blobfiltPyramid=load([folderData,'\blobfiltPyramid.mat']);
blobfiltPyramid=blobfiltPyramid.blobfiltPyramid;

% Load 2D data
blobfiltPyramid2D_Raw=load([folderData,'\blobfilt2D.mat']);
blobfiltPyramid2D_Raw=blobfiltPyramid2D_Raw.blobfilt2D;

% Load kmeans3D
kmeans3D=load([folderData,'\kmeans3D.mat']);
kmeans3D=kmeans3D.kmeans3D;

pyramidalNeurons=kmeans3D.PyramidalNeurons;
blobfiltPyramid2D= blobfiltPyramid2D_Raw(pyramidalNeurons, :); % Table with only artefacts
pyraRadius2D =blobfiltPyramid2D.SterelogyDistance;

window =load([folderData,'\window.mat']);
window=window.window;
windowRes = round(window*voxelRes/1000); % convert window to µm dimensions

end

function [cellTable]=plotEstimationFnc(folder,windowRes,name,voxelRes,z,blobfiltPyramid,pyraRadius2D)
% Plot and save data
% Sintax:
%     plotEstimationFnc(folder,blobfiltPyramid,windowRes,name,voxelRes,z)
% Inputs:
%     folder,           Folder name for saved images
%     blobfiltPyramid,  Table with information of detected pyramidal cells
%     windowRes,        Defined sampling window adjusted with the voxel
%                       resolution
%     name,             name of subject/folder that will be analysed
%     voxelRes,         Voxel resolution of the analysed data
%     z,                Z-height between each section

%  Create tables
cellTable.cell_vol =blobfiltPyramid.Volume*voxelRes/1000*voxelRes/1000*z/1000;% Estimate Volume based on number of voxels and resolution
cellTable.cell_vol_mean= round(mean(cellTable.cell_vol,'omitnan'));
cellTable.cell_vol_std= round(std(cellTable.cell_vol,'omitnan'));
cellTable.cell_vol_diameter= nthroot(3*cellTable.cell_vol_mean/(4*pi),3)*2;
cellTable.cell_vol_radius= nthroot(3*cellTable.cell_vol_mean/(4*pi),3);
cellTable.cell_volLog = (log(cellTable.cell_vol));
cellTable.cell_vol_log_mean =round(exp(mean(cellTable.cell_volLog)));
cellTable.cell_vol_log_std = round(exp(std(cellTable.cell_volLog)));

%%%%%%%%%%%%%% Log normal distribution
% distributionFitter(cellTable.cell_vol)
pd=fitdist(cellTable.cell_vol,'Lognormal'); % Fit data to distribution
cellTable.cell_vol_logFit_mean = round(mean(pd));
cellTable.cell_vol_logFit_std = round(var(pd));
cellTable.cell_vol_logFit_mode = round(exp(pd.mu-pd.sigma^2));

cellTable.cell_Sphericity = blobfiltPyramid.Sphericity;
cellTable.cell_Sphericity_mean =round(mean(cellTable.cell_Sphericity,'omitnan')*100)/100;
cellTable.cell_Sphericity_std = round(std(cellTable.cell_Sphericity,'omitnan')*100)/100;
cellTable.cell_SphericityLog = (log(cellTable.cell_Sphericity));
cellTable.cell_Sphericity_log_mean =round(exp(mean(cellTable.cell_SphericityLog,'omitnan'))*100)/100;
cellTable.cell_Sphericity_log_std = (exp((mean(cellTable.cell_SphericityLog))));

% Diameter
cellTable.cell_diameter = pyraRadius2D*2;
cellTable.cell_diameter_mean =round(mean(cellTable.cell_diameter)*100)/100;
cellTable.cell_diameter_std = round(std(cellTable.cell_diameter)*100)/100;
cellTable.cell_diameterLog = (log(cellTable.cell_diameter));
cellTable.cell_diameter_log_mean =round(exp(mean(cellTable.cell_diameterLog))*100)/100;
cellTable.cell_diameter_log_std = round(exp(std(cellTable.cell_diameterLog))*100)/100;
cellTable.cell_diameter_Vol =round(4/3*pi*(mean(pyraRadius2D)).^3); % Estimate volume
cellTable.cell_radius = mean(pyraRadius2D,'omitnan');

% Orientation
cellTable.cell_orientation = (blobfiltPyramid.thetaFeret);
cellTable.cell_orientation_mean =round(mean(cellTable.cell_orientation,'omitnan'));
cellTable.cell_orientation_std = round(std(cellTable.cell_orientation,'omitnan'));
cellTable.cell_orientationLog = log(cellTable.cell_orientation);
cellTable.cell_orientationLog(isinf(cellTable.cell_orientationLog)|isnan(cellTable.cell_orientationLog)) = [];
cellTable.cell_orientation_log_mean =round(exp(mean(cellTable.cell_orientationLog ,'omitnan')));
cellTable.cell_orientation_log_std = round(exp(std(cellTable.cell_orientationLog ,'omitnan')));

% homogeneity in x, y and z-axis of centroids
cellTable.cell_x = blobfiltPyramid.Centroid(:,1);
cellTable.cell_y = blobfiltPyramid.Centroid(:,2);
cellTable.cell_z = blobfiltPyramid.Centroid(:,3);
cellTable.density = round(size(blobfiltPyramid,1)/(windowRes(1)*windowRes(2)*windowRes(3))*10^9); % Density for mm^3
%%
%%%%%%%%%%%%%% Figures
close all;
N = 50;
figure(1);
edges = linspace(0,max(cellTable.cell_vol),N); % Create 20 bins.
histogram(cellTable.cell_vol, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
xlabel('Volume (µm^3)','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',12);
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
set(gca, 'XLim', [0, get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['VolumePlot_',name]),'m')
saveas(gcf,fullfile(folder, ['VolumePlot_',name]),'pdf')

figure(2)
edges = linspace(min(cellTable.cell_volLog),max(cellTable.cell_volLog),N); % Create 20 bins.
histogram(cellTable.cell_volLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',4,'FontWeight','bold')
xticks((linspace(min(cellTable.cell_volLog),max(cellTable.cell_volLog),7)));
xlabel('Log Volume (µm^3)','FontSize',4,'FontWeight','bold')
set(gca,'FontSize',12);
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlim([0 max(((cellTable.cell_volLog)))]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round(exp(xtix)))
set(gca, 'XLim', [round(min(cellTable.cell_volLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['VolumePlotLog_',name]),'m')
saveas(gcf,fullfile(folder, ['VolumePlotLog_',name]),'pdf')

figure(3)
edges = linspace(min(cellTable.cell_Sphericity),max(cellTable.cell_Sphericity),N); % Create 20 bins.
histogram(cellTable.cell_Sphericity, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',12,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('Sphericity','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
saveas(gcf,fullfile(folder, ['SphericityPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['SphericityPlot_',name]),'pdf')

blobfiltPyramid.Sphericity(blobfiltPyramid.Sphericity>0.37)=blobfiltPyramid.Sphericity(blobfiltPyramid.Sphericity>0.37)-0.13;

figure(4)
edges = linspace(min(cellTable.cell_SphericityLog),max(cellTable.cell_SphericityLog),N); % Create 20 bins.
histogram(cellTable.cell_SphericityLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('Log Sphericity','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xlim([min(((cellTable.cell_SphericityLog))) max(((cellTable.cell_SphericityLog)))]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round(((exp(xtix))*100))/100)
set(gca, 'XLim', [round(min(cellTable.cell_SphericityLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['SphericityPlotLog_',name]),'m')
saveas(gcf,fullfile(folder, ['SphericityPlotLog_',name]),'pdf')

figure(5)
edges = linspace(0,max(cellTable.cell_orientation),N); % Create 20 bins.
histogram(cellTable.cell_orientation, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1);
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlim([0 max(cellTable.cell_orientation)]);
xlabel(['Orientation (',char(176),')'],'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
% ax = hINT.Parent;   % Important
% set(ax, 'XTick', 0:10:90)
saveas(gcf,fullfile(folder, ['OrientationPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['OrientationPlot_',name]),'pdf')

figure(6)
edges = linspace(0,max(cellTable.cell_orientationLog),N); % Create 20 bins.
histogram(cellTable.cell_orientationLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1);
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlim([0 max(((cellTable.cell_orientationLog)))]);
xlabel(['Log Orientation (',char(176),')'],'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((exp(xtix))))
set(gca, 'XLim', [round(min(cellTable.cell_orientationLog),3),get(gca, 'XLim') *[0; 1]])
saveas(gcf,fullfile(folder, ['OrientationLogPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['OrientationLogPlot_',name]),'pdf')

figure(7)
edges = linspace(min(cellTable.cell_diameter),max(cellTable.cell_diameter),N); % Create 20 bins.
histogram(cellTable.cell_diameter, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
% xlim([0 max(cellTable.cell_orientation)]);
xlabel('Diameter (µm)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((xtix)))
set(gca, 'XLim', [round(min(cellTable.cell_diameterLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['DiameterPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['DiamterPlot_',name]),'pdf')

figure(8)
edges = linspace(min(cellTable.cell_diameterLog),max(cellTable.cell_diameterLog),N); % Create 20 bins.
histogram(cellTable.cell_diameterLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
% xlim([0 max(cellTable.cell_orientation)]);
xlabel('Log Diameter (µm)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xlim([0 max(((cellTable.cell_diameterLog)))]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((exp(xtix))))
set(gca, 'XLim', [round(min(cellTable.cell_diameterLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['DiameternLogPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['DiamterLogPlot_',name]),'pdf')

save(fullfile(folder, ['cellTable_',name,'.mat']), 'cellTable')
end
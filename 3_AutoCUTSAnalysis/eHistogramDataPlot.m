% Code to segmentation of sections.
% Code written and posted by Nick Yin Larsen, November 2020.
%------------------------------------------------------------------------------------------------
% Startup code.
% Measured properties of the pyramidal neurons
clear all;
voxelRes = 272; % nm
z = 800; % nm
name = 'subject1';
folder1 ='AutuCUTS_Pipeline';
folder2 ='Example_4_Prediction';
folder3 ='Example_UNetDense_stack_layer3_results';

[folder,blobfiltPyramid,windowRes]=readFilesFnc(folder1,folder2,folder3,name,voxelRes);
plotEstimationFnc(folder,blobfiltPyramid,windowRes,name,voxelRes,z)

%% %%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%
function [folder,blobfiltPyramid,windowRes]=readFilesFnc(folder1,folder2,folder3,name,voxelRes)
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
folderSave =[name,'_histogram_results'];
folder = [savePath,'\',folderSave];
folderNameSave=folder;
if ~exist(folderNameSave, 'dir')
    mkdir(folderNameSave);
end

% Load data
folderData = [savePath,'\',folder3];
load([folderData,'\blobfiltPyramid.mat']);
load([folderData,'\window.mat']);
windowRes = round(window*voxelRes/1000); % convert window to µm dimensions

end

function plotEstimationFnc(folder,blobfiltPyramid,windowRes,name,voxelRes,z)
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
cellTable.cell_volLog = log10(cellTable.cell_vol);
cellTable.cell_vol_median =median(cellTable.cell_vol);
VolQ25 = quantile(cellTable.cell_vol,0.25);
VolQ75 = quantile(cellTable.cell_vol,0.75);
cellTable.cell_vol_spread =  (VolQ75-VolQ25)/2;
cellTable.cell_vol_log_mean =mean(10.^cellTable.cell_volLog);
cellTable.cell_vol_log_std = std(10.^cellTable.cell_volLog);
% http://davidmlane.com/hyperstat/A48607.html

cellTable.cell_Sphericity = blobfiltPyramid.Sphericity;
cellTable.cell_SphericityLog = log10(cellTable.cell_Sphericity);
cellTable.cell_Sphericity_mean =mean(cellTable.cell_Sphericity);
cellTable.cell_Sphericity_std = std(cellTable.cell_Sphericity);
cellTable.cell_Sphericity_log_mean =mean(10.^cellTable.cell_SphericityLog);
cellTable.cell_Sphericity_log_std = std(10.^cellTable.cell_SphericityLog);

% Diameter
cellTable.cell_diameter = blobfiltPyramid.Radius*2;
cellTable.cell_diameter_mean =mean(cellTable.cell_diameter);
cellTable.cell_diameter_std = std(cellTable.cell_diameter);
cellTable.cell_diameterLog = log10(cellTable.cell_diameter);
cellTable.cell_diameter_log_mean =mean(10.^cellTable.cell_diameterLog);
cellTable.cell_diameter_log_std = std(10.^cellTable.cell_diameterLog);

% Orientation
filterOri = blobfiltPyramid.theta>90;        %% Remove angles above 90 degrees
cellTable.cell_orientation = (blobfiltPyramid.theta(~filterOri));
cellTable.cell_orientationLog = log10(cellTable.cell_orientation);
cellTable.cell_orientation_median =median(cellTable.cell_orientation);
oriQ25 = quantile(cellTable.cell_orientation,0.25);
oriQ75 = quantile(cellTable.cell_orientation,0.75);
cellTable.cell_orientation_spread =  (oriQ75-oriQ25)/2;
cellTable.cell_orientation_mean =mean(cellTable.cell_orientation);
cellTable.cell_orientation_std = std(cellTable.cell_orientation);
cellTable.cell_orientation_log_mean =mean(10.^cellTable.cell_orientationLog);
cellTable.cell_orientation_log_std = std(10.^cellTable.cell_orientationLog);


cellTable.cell_x = blobfiltPyramid.Centroid(:,1);
cellTable.cell_y = blobfiltPyramid.Centroid(:,2);
cellTable.cell_z = blobfiltPyramid.Centroid(:,3);
cellTable.density = round(size(blobfiltPyramid,1)/(windowRes(1)*windowRes(2)*windowRes(3))*10^9); % Density for mm^3
save(fullfile(folder, ['cellTable_',name,'.mat']), 'cellTable')
%%
%%%%%%%%%%%%%% Figures
close all;
N = 25;
figure;
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
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((10.^xtix)))
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


figure(4)
edges = linspace(min(cellTable.cell_SphericityLog),max(cellTable.cell_SphericityLog),N); % Create 20 bins.
histogram(cellTable.cell_SphericityLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('Log Sphericity','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((10.^xtix)*100)/100)
set(gca, 'XLim', [round(min(cellTable.cell_SphericityLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['SphericityPlotLog_',name]),'m')
saveas(gcf,fullfile(folder, ['SphericityPlotLog_',name]),'pdf')

figure(5)
edges = linspace(0,max(cellTable.cell_orientation),N); % Create 20 bins.
histogram(cellTable.cell_orientation, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlim([0 max(cellTable.cell_orientation)]);
xlabel(['Orientation (',char(176),')'],'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
saveas(gcf,fullfile(folder, ['OrientationPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['OrientationPlot_',name]),'pdf')

figure(6)
edges = linspace(min(cellTable.cell_orientationLog),max(cellTable.cell_orientationLog),N); % Create 20 bins.
histogram(cellTable.cell_orientationLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel(['Log Orientation (',char(176),')'],'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((10.^xtix)))
set(gca, 'XLim', [round(min(cellTable.cell_orientationLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['OrientationLogPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['OrientationLogPlot_',name]),'pdf')

figure(7)
edges = linspace(min(cellTable.cell_diameterLog),max(cellTable.cell_diameterLog),N); % Create 20 bins.
histogram(cellTable.cell_diameterLog, 'Normalization','probability', 'BinEdges',edges,'FaceColor','none', 'LineWidth', 1)
ylabel('% of Neurons','FontSize',14,'FontWeight','bold')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
% xlim([0 max(cellTable.cell_orientation)]);
xlabel(['Log Diamter (µm)'],'FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',round((10.^xtix)))
set(gca, 'XLim', [round(min(cellTable.cell_diameterLog),3), get(gca, 'XLim') * [0; 1]])
saveas(gcf,fullfile(folder, ['DiameternLogPlot_',name]),'m')
saveas(gcf,fullfile(folder, ['DiamterLogPlot_',name]),'pdf')
end
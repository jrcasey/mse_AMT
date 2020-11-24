%% Draw AMT13 cruise transect on a map

% uses m_map: https://www.eoas.ubc.ca/~rich/map.html

%% Load data
load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/output/FullSolution_L2.mat');
FullSolution = FullSolution_L2;

%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;

%% 
fig = figure('position',[500 500 500 1000])

m_proj('lambert','long',[-70 0],'lat',[-50 60]);
[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -500 -200 0 ],'edgecolor','none');
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
h2=m_line(CruiseData.Lon(Gridding.stationsVec2), CruiseData.Lat(Gridding.stationsVec2),'marker','o','color','r','linewi',2,...
          'linest','none','markersize',8,'markerfacecolor','w');
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
%legend(h2(1),'AMT-13 Stations','location','southwest');


colormap(m_colmap('blues'));  
caxis([-7000 000]);

%[ax,h]=m_contfbar([.55 .75],.8,CS,CH,'endpiece','no','axfrac',.05);
%title(ax,'meters')

set(gcf,'color','w');  % otherwise 'print' turns lakes black

saveas(fig,'/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/CruiseVariables/Transect.eps','epsc')
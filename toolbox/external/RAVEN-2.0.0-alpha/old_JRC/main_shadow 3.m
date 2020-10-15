%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For John CASY taken on 17/09/2014
%% Edited on 18/09/2014
%% Partho sen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

disp('load model draftMED4...')
load('model.mat');

%model = importModel('draftMED4_v2.13.xml');
%%exportToExcelFormat(model,'draftMED4.xls');

%[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

%% GET BIOMASS REACTIONS

disp('Get Biomass reactions...')
ind = find(strcmp(model.rxnNames,'Biomass formation'));

%RXNS = model.rxnNames(ind);

%% SET OBJECTIVE FUNCTION

disp('Set Objective function to 1...')
model.c(ind) = 1; %% set objective function as 'biomass of formation'


%%% SIMULATE RAVEN %% file written

% try applying a stressor. Nitrogen:
%model = setParam(model,'lb','AmmoniaTRANS',-1000); % 0.6 is optimal
%model = setParam(model,'eq','UreaTRANS',0);
%model = setParam(model,'eq','CyanateTRANS',0);

% try applying a second stressor. Light:
%model = setParam(model,'eq','LightTRANS',1); % 52 is optimal

% and phosphate
%model = setParam(model,'lb','OrthophosphateTRANS',-1000); % -0.03 is optimal

disp('Simulate model MOSEK(RAVEN)...')

[Soln_RAV,hsSolOut,res] = SolveLP2(model,'max');

disp('Plot Shadow price...')

[sneg_idx,spos_idx] = plot_shadow_RAV(model,Soln_RAV,1001,'false'); %% 'TRUE' writes a File in the same directory


%%% SIMULATE COBRA %% file written

% disp('Simulate model GLPX(COBRA)...')
% 
% [Soln_COB] = optimizeCbModel(model,'max');
% 
% disp('Plot Shadow price...')
% 
% [sneg_idx1,spos_idx1] = plot_shadow_COB(model,Soln_COB,1002,'false'); %% 'TRUE' writes a File in the same directory
% 


return

%%%% The END %%%%%%








%%% Changing bounds

%  model.b(sneg_idx) = ones(1,length(sneg_idx)).*(0); %%% on pertubing b no optimal soluation found
%  Solnyy = solveLP(model,'max');
% 
% [sneg_idx,spos_idx] = plot_shadow(model,Solnyy,1002);
% 










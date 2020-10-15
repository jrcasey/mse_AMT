%% Generate model from arrayData
function INITmodel = getINITmodelFromArrayDataInServer(filename,filename_index,scoring_scale)
% getINITmodelFromArrayDataInServer(filename_index,scoring_scale)
%
% This function generates an INIT model from ArrayData in a file. refModel, 
% and tasks must be specified within this script. scoring_scale defines
% either a stringent (1) or a relaxed (2) score.

%%% Initialize path
ravenPath  = genpath('/c3se/users/gatto/Glenn/RAVEN/');
mosek6path = genpath('/c3se/users/gatto/Glenn/mosek/6');
addpath(ravenPath);
addpath(mosek6path);

%%% Load Data
load iTO977.mat                     %Reference model
load tasks_iTO997_20140423.mat      %Tasks
if exist('hotStart.mat','file')==2
    load hotStart.mat                   %MILP Hot Start
    hotStart              = res.sol.int.xx;
else
    hotStart = [];
end
if exist('essentialRxnsForTasks_iTO997.mat','file')==2
    load essentialRxnsForTasks_iTO997.mat         %Essential rxns array
else
    [~,essential] = checkTasks(model,[],true,false,true,taskList);
    essentialRxnsForTasks = model.rxns(any(essential,2));
    save('essentialRxnsForTasks_iTO997','essentialRxnsForTasks')
end

%%% Files
refModel              = model;
tissue                = 'S. cerevisiae';
celltype              = filename{filename_index}(12:end);
hpaData               = [];
arrayData             = parseArray(filename{filename_index},tissue,celltype,{'RNAseq'},{'Supportive'});
metabolomicsData      = [];

%%% Settings
allowExcretion      = true;  %Enforce presence of reactions that may be regarded as dead-ends just because a complete biomass eqn is missing
printReport         = true;
hpaLevelScores.names   = {};
hpaLevelScores.scores  = [];
if scoring_scale == 1
    arrayDataScores.names  = {'Very high' 'High' 'Medium' 'Low' 'Very low' 'None'};
    arrayDataScores.scores = [20 10 2 -2 -10 -20];
elseif scoring_scale == 2
    arrayDataScores.names  = {'Very high' 'High' 'Medium' 'Low' 'Very low' 'None'};
    arrayDataScores.scores = [20 10 2 0 -2 -20];
end
params.maxTime      = 96*60;
params.printReport  = true;
params.relGap       = 0.001;

%%% Run tINIT
fprintf('Running INIT...\n')
[INITmodel metProduction deletedDeadEndRxns deletedRxnsInINIT] = ...
    getINITModelInServer(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, ...
    printReport,params, allowExcretion, hpaLevelScores,arrayDataScores,essentialRxnsForTasks,hotStart);
fprintf('\tINIT successfully terminated. Saving output...\n')
outputfile = strcat('INITresults_',filename{filename_index}(12:end));
    save(outputfile,'INITmodel','deletedDeadEndRxns','deletedRxnsInINIT')
fprintf('Saved.\n')
end

%This is to get only INIT using mosek 6
function [initModel metProduction deletedDeadEndRxns deletedRxnsInINIT]=getINITModelInServer(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, printReport, params, allowExcretion, hpaLevelScores,arrayDataScores,essentialRxnsForTasks,hotStart)
% getINITModelInServer
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%
%   refModel            a model structure. The model should be in the
%                       closed form (no exchange reactions open). Import
%                       using import(filename,false). If the model is not
%                       loaded using importModel, it might be that there
%                       is no "unconstrained" field. In that case,
%                       manually add the field like:
%                       model.unconstrained=false(numel(model.mets),1);
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this 
%                       tissue (opt, default is to use the best values
%                       among all the cell types for the tissue. Use [] if
%                       you want to supply more arguments)
%   hpaData             HPA data structure from parseHPA (opt if arrayData is 
%                       supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not be
%                       unique, as there can be multiple cell types per tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%   metabolomicsData    cell array with metabolite names that the model
%                       should produce (opt, default [])
%   taskFile            a task list in Excel format. See parseTaskList for
%                       details (opt, default [])
%   useScoresForTasks   true if the calculated reaction scored should be used as
%                       weights in the fitting to tasks (opt, default true)
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   taskStructure       task structure as from parseTaskList. Can be used
%                       as an alternative way to define tasks when Excel
%                       sheets are not suitable. Overrides taskFile (opt,
%                       default [])
%   params              parameter structure as used by getMILPParams. This is
%                       for the INIT algorithm. For the the MILP problems
%                       solved to fit tasks, see paramsFT (opt, default [])
%   paramsFT            parameter structure as used by getMILPParams. This is
%                       for the fitTasksInServer step. For the INIT algorithm, see
%                       params (opt, default [])
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                     This results in fewer reactions being considered
%                     dead-ends, but all reactions in the resulting model may
%                     not be able to carry flux. If this is "false" then the
%                     equality constraints are taken from model.b. If the
%                     input model lacks exchange reactions then this should
%                     probably be "true", or a large proportion of the model
%                     would be excluded for connectivity reasons
%                     (opt, default false)
%
%   model                   the resulting model structure
%   metProduction           array that indicates which of the
%                           metabolites in metabolomicsData that could be
%                           produced. Note that this is before the
%                           gap-filling process to enable defined tasks. To
%                           see which metabolites that can be produced in
%                           the final model, use canProduce. 
%                           -2: metabolite name not found in model
%                           -1: metabolite found, but it could not be produced
%                           1: metabolite could be produced
%   essentialRxnsForTasks   cell array of the reactions which were
%                           essential to perform the tasks
%   addedRxnsForTasks       cell array of the reactions which were added in
%                           order to perform the tasks
%   deletedDeadEndRxns      cell array of reactions deleted because they
%                           could not carry flux (INIT requires a
%                           functional input model)
%   deletedRxnsInINIT       cell array of the reactions which were deleted by 
%                           the INIT algorithm
%   taskReport              structure with the results for each task
%   	id                  cell array with the id of the task
%       description         cell array with the description of the task
%       ok                  boolean array with true if the task was successful
%       essential           cell array with cell arrays of essential
%                           reactions for the task
%       gapfill             cell array of cell arrays of reactions included
%                           in the gap-filling for the task
%
%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreModel, runINIT and fitTasksInServer individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasksInServer rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [model metProduction essentialRxnsForTasks addedRxnsForTasks...
%               deletedDeadEndRxns deletedRxnsInINIT taskReport]=...
%               getINITModel(refModel, tissue, celltype, hpaData, arrayData,...
%               metabolomicsData, taskFile, useScoresForTasks, printReport,...
%               taskStructure, params, paramsFT)
%
%   Rasmus Agren, 2013-05-16
%

if nargin<3
    celltype=[];
end
if nargin<4
    hpaData=[];
end
if nargin<5
    arrayData=[];
end
if nargin<6
    metabolomicsData=[];
end
if nargin<7
    printReport=true;
end
if nargin<8
    params=[];
end
if nargin<9
    allowExcretion=false;
end
if nargin<10
    hpaLevelScores.names={'High' 'Medium' 'Low' 'None' 'Strong' 'Moderate' 'Weak' 'Negative'};
    hpaLevelScores.scores=[20 15 10 -8 20 15 10 -8];
end
if nargin<11
    arrayDataScores.names={'Very high' 'High' 'Medium' 'Low' 'Very low' 'None'};
    arrayDataScores.scores=[20 10 2 -2 -10 -20];
end
if nargin < 12
    essentialRxnsForTasks = {};
end
if nargin < 13
    hotStart = [];
end

%Check that the model is in the closed form
if ~isfield(refModel,'unconstrained')
   throw(MException('','Exchange metabolites should normally not be removed from the model when using getINITModel. Use importModel(file,false) to import a model with exchange metabolites remaining (see the documentation for details)')); 
end

if printReport==true
    if any(celltype)
        fprintf(['***Generating model for: ' tissue ' - ' celltype '\n']);
    else
        fprintf(['***Generating model for: ' tissue '\n']);
    end
    if ~isempty(hpaData)
        fprintf('-Using HPA data\n');
    end
    if ~isempty(arrayData)
        fprintf('-Using array data\n');
    end
    if ~isempty(metabolomicsData)
        fprintf('-Using metabolomics data\n');
    end
    fprintf('\n');
    
    printScores(refModel,'Reference model statistics',hpaData,arrayData,tissue,celltype);
end


%Remove dead-end reactions to speed up the optimization and to
%differentiate between reactions removed by INIT and those that are
%dead-end
[crap, deletedDeadEndRxns]=simplifyModel(refModel,true,false,true,true,true);
cModel=removeRxns(refModel,deletedDeadEndRxns,false,true);

%Store the connected model like this to keep track of stuff
if printReport==true
    printScores(cModel,'Pruned model statistics',hpaData,arrayData,tissue,celltype);
end

%Score the connected model
noGeneScore         = 0;
multipleGeneScoring = 'best';
multipleCellScoring = 'best';
[rxnScores geneScores]=scoreModelwArray(cModel,hpaData,arrayData,tissue,celltype,noGeneScore,multipleGeneScoring,multipleCellScoring,hpaLevelScores,arrayDataScores);

%Run the INIT algorithm. The exchange reactions that are used in the final
%reactions will be open, which doesn't fit with the last step. Therefore
%delete reactions from the original model instead of taking the output.
%The default implementation does not constrain reversible reactions to only
%carry flux in one direction.
%Runs without the constraints on reversibility and with all output allowed.
%This is to reduce the complexity of the problem.

[crap deletedRxnsInINIT metProduction]=runINITinServer(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,allowExcretion,false,params,hotStart);
initModel=removeRxns(cModel,deletedRxnsInINIT,true,true);
if printReport==true
    printScores(initModel,'INIT model statistics',hpaData,arrayData,tissue,celltype);
    printScores(removeRxns(cModel,setdiff(cModel.rxns,deletedRxnsInINIT),true,true),'Reactions deleted by INIT',hpaData,arrayData,tissue,celltype);
end
%The full model has exchange reactions in it. fitTasksInServer calls on fillGaps,
%which automatically removes exchange metabolites (because it assumes that
%the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are delete from the
%reference model here.
initModel.id='INITModel';
end

%This is for printing a summary of a model
function [rxnS geneS]=printScores(model,name,hpaData,arrayData,tissue,celltype)
    [a b]=scoreModelwArray(model,hpaData,arrayData,tissue,celltype);
    rxnS=mean(a);
    geneS=mean(b(~isinf(b)));
    fprintf([name ':\n']);
    fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
    fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
    fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
    fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end

%This is to use array structure to score the rxns
function [rxnScores geneScores hpaScores arrayScores]=scoreModelwArray(model,hpaData,arrayData,tissue,celltype,noGeneScore,multipleGeneScoring,multipleCellScoring,hpaLevelScores,arrayDataScores)
% scoreModelwArray
%   Scores the reactions and genes in a model based on expression data
%   from HPA and/or gene arrays
%
%   model               a model structure
%   hpaData             HPA data structure from parseHPA (opt if arrayData is 
%                       supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not be
%                       unique, as there can be multiple cell types per tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this 
%                       tissue (opt, default is to use the best values
%                       among all the cell types for the tissue. Use [] if
%                       you want to supply more arguments)
%   noGeneScore         score for reactions without genes (opt, default -2)
%   multipleGeneScoring determines how scores are calculated for reactions
%                       with several genes ('best' or 'average')
%                       (opt, default 'best')
%   multipleCellScoring determines how scores are calculated when several
%                       cell types are used ('best' or 'average')
%                       (opt, default 'best')
%   hpaLevelScores      structure with numerical scores for the expression
%                       level categories from HPA. The structure should have a
%                       "names" and a "scores" field (opt, see code for
%                       default scores)
%
%   rxnScores       scores for each of the reactions in model
%   geneScores      scores for each of the genes in model. Genes which are
%                   not in the dataset(s) have -Inf as scores
%   hpaScores       scores for each of the genes in model if only taking hpaData
%                   into account. Genes which are not in the dataset(s) 
%                   have -Inf as scores
%   arrayScores     scores for each of the genes in model if only taking arrayData
%                   into account. Genes which are not in the dataset(s) 
%                   have -Inf as scores
%       
%   Usage: [rxnScores geneScores hpaScores arrayScores]=scoreModel(model,...
%               hpaData,arrayData,tissue,celltype,noGeneScore,multipleGeneScoring,...
%               multipleCellScoring,hpaLevelScores)
%
%   Rasmus Agren, 2012-09-11
%

if nargin<3
    arrayData=[];
end
if nargin<5
    celltype=[];
end
if nargin<6
    noGeneScore=0;
end
if nargin<7
    multipleGeneScoring='best';
end
if nargin<8
    multipleCellScoring='best';
end
if nargin<9
    %The first four are for APE, the other ones for staining
    hpaLevelScores.names={'High' 'Medium' 'Low' 'None' 'Strong' 'Moderate' 'Weak' 'Negative'};
    hpaLevelScores.scores=[20 15 10 -8 20 15 10 -8];
end
if nargin<10
    arrayDataScores.names={'Very high' 'High' 'Medium' 'Low' 'Very low' 'None'};
    arrayDataScores.scores=[20 10 2 -2 -10 -20];
end

if isempty(hpaData) && isempty(arrayData)
    throw(MException('','Must supply hpaData, arrayData or both'));    
end
if ~strcmpi(multipleGeneScoring,'best') && ~strcmpi(multipleGeneScoring,'average')
    throw(MException('','Valid options for multipleGeneScoring are "best" or "average"'));    
end
if ~strcmpi(multipleCellScoring,'best') && ~strcmpi(multipleCellScoring,'average')
    throw(MException('','Valid options for multipleCellScoring are "best" or "average"'));    
end

%This is so that the code can ignore which combination of input data that is
%used
if isempty(arrayData)
    usingArray = false;
    arrayData.genes={};
    arrayData.tissues={};
    arrayData.celltypes={};
    arrayData.levels={};
    arrayData.types={};
    arrayData.reliabilities={};
    arrayData.gene2Level=[];
    arrayData.gene2Type=[];
    arrayData.gene2Reliability=[];
else
    usingArray = true;
end
if isempty(hpaData)
    usingHpa = false;
    hpaData.genes={};
    hpaData.tissues={};
    hpaData.celltypes={};
    hpaData.levels={};
    hpaData.types={};
    hpaData.reliabilities={};
    hpaData.gene2Level=[];
    hpaData.gene2Type=[];
    hpaData.gene2Reliability=[];
else
    usingHpa = true;
end

%Check that the tissue exists
if usingArray
    if usingHpa
        if ~ismember(upper(tissue),upper(hpaData.tissues)) && ~ismember(upper(tissue),upper(arrayData.tissues))
            throw(MException('','The tissue name does not match')); 
        end
    else
        if ~ismember(upper(tissue),upper(arrayData.tissues))
            throw(MException('','The tissue name does not match'));  
        end
    end
else
    if usingHpa
        if ~ismember(upper(tissue),upper(hpaData.tissues))
            throw(MException('','The tissue name does not match')); 
        end
    end    
end

if any(celltype)
    %Check that both data types has cell type defined if that is to be used
    if ~isfield(hpaData,'celltypes') || ~isfield(arrayData,'celltypes')
        throw(MException('','Both hpaData and arrayData must contain cell type information if cell type is to be used'));   
    end
    if usingArray
        if usingHpa
            if ~ismember(upper(celltype),upper(hpaData.celltypes)) && ~ismember(upper(celltype),upper(arrayData.celltypes))
                throw(MException('','The cell type name does not match'));   
            end
        else
            if ~ismember(upper(celltype),upper(arrayData.celltypes))
                throw(MException('','The cell type name does not match'));   
            end
        end
    else
        if usingHpa
            if ~ismember(upper(celltype),upper(hpaData.celltypes))
                throw(MException('','The cell type name does not match'));   
            end
        end    
    end
end

%Some preprocessing of the structures to speed up a little
%Remove all tissues that are not the correct one
J=~strcmpi(hpaData.tissues,tissue);

%If cell type is supplied, then only keep that cell type
if any(celltype)
    J=J | ~strcmpi(hpaData.celltypes,celltype);
end

hpaData.tissues(J)=[];
if isfield(hpaData,'celltypes')
    hpaData.celltypes(J)=[];
end
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(:,J)=[];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(:,J)=[];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(:,J)=[];
end

%Remove all genes from the structures that are not in model or that aren't
%measured in the tissue
if ~isempty(hpaData.genes) %This should not be necessary, but the summation is a 0x1 matrix and the other is []
    I=~ismember(hpaData.genes,model.genes) | sum(hpaData.gene2Level,2)==0;
else
    I=[];
end
hpaData.genes(I)=[];
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(I,:)=[];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(I,:)=[];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(I,:)=[];
end

%Some preprocessing of the structures to speed up a little
%Remove all tissues that are not the correct one
J=~strcmpi(arrayData.tissues,tissue);

%If cell type is supplied, then only keep that cell type
if any(celltype)
    J=J | ~strcmpi(arrayData.celltypes,celltype);
end

arrayData.tissues(J)=[];
if isfield(arrayData,'celltypes')
    arrayData.celltypes(J)=[];
end
if isfield(arrayData,'gene2Level')
    arrayData.gene2Level(:,J)=[];
end
if isfield(arrayData,'gene2Type')
    arrayData.gene2Type(:,J)=[];
end
if isfield(arrayData,'gene2Reliability')
    arrayData.gene2Reliability(:,J)=[];
end

%Remove all genes from the structures that are not in model or that aren't
%measured in the tissue
if ~isempty(arrayData.genes) %This should not be necessary, but the summation is a 0x1 matrix and the other is []
    I=~ismember(arrayData.genes,model.genes) | sum(arrayData.gene2Level,2)==0;
else
    I=[];
end
arrayData.genes(I)=[];
if isfield(arrayData,'gene2Level')
    arrayData.gene2Level(I,:)=[];
end
if isfield(arrayData,'gene2Type')
    arrayData.gene2Type(I,:)=[];
end
if isfield(arrayData,'gene2Reliability')
    arrayData.gene2Reliability(I,:)=[];
end

%Map the array levels to scores
[I J]=ismember(upper(arrayData.levels),upper(arrayDataScores.names));
if ~all(I)
    throw(MException('','There are expression level categories that do not match to hpaLevelScores'));  
end
[K L M]=find(arrayData.gene2Level);
scores=arrayDataScores.scores(J);
if strcmpi(multipleCellScoring,'best')
    aScores=max(sparse(K,L,scores(M),numel(arrayData.genes),numel(arrayData.tissues)),[],2);
else
    aScores=mean(sparse(K,L,scores(M),numel(arrayData.genes),numel(arrayData.tissues)),2);
end
% 
% %Calculate the scores for the arrayData. These scores are
% %calculated for each genes from its fold change between the tissue/celltype(s) in
% %question and all other celltypes. This is a lower quality data than
% %protein abundance, since a gene that is equally highly expressed in all cell types
% %will have a score of 0.0. These scores are therefore only used for genes
% %for which there is no HPA data available. The fold changes are transformed
% %as min(log(x),10) for x>1 and max(log(x),-5) for x<1 in order to have
% %negative scores for lower expressed genes and to scale the scrores to have
% %somewhat lower weights than the HPA scores
% tempArrayLevels=arrayData.levels;
% tempArrayLevels(isnan(tempArrayLevels))=0;
% average=sum(tempArrayLevels,2)./sum(~isnan(arrayData.levels),2);
% if strcmpi(multipleCellScoring,'best')
%     current=max(tempArrayLevels(:,I),[],2);
% else
%     current=sum(tempArrayLevels(:,I),2)./sum(~isnan(arrayData.levels(:,I)),2);
% end
% if size(current,2) > 1
%     if any(current)
%         aScores=5*log(current./average);
%     else
%         aScores=[];
%     end
%     aScores(aScores>0)=min(aScores(aScores>0),10);
%     aScores(aScores<0)=max(aScores(aScores<0),-5);
% else
%     aScores= 20*(current>1) -20*(current<-1);
% end

%Map the HPA levels to scores
[I J]=ismember(upper(hpaData.levels),upper(hpaLevelScores.names));
if ~all(I)
    throw(MException('','There are expression level categories that do not match to hpaLevelScores'));  
end
[K L M]=find(hpaData.gene2Level);
scores=hpaLevelScores.scores(J);
if strcmpi(multipleCellScoring,'best')
    hScores=max(sparse(K,L,scores(M),numel(hpaData.genes),numel(hpaData.tissues)),[],2);
else
    hScores=mean(sparse(K,L,scores(M),numel(hpaData.genes),numel(hpaData.tissues)),2);
end

%Get the scores for the genes, only use HPA if available
geneScores=inf(numel(model.genes),1)*-1;
hpaScores=geneScores;
arrayScores=geneScores;

[I J]=ismember(model.genes,hpaData.genes);
hpaScores(I)=hScores(J(I));
geneScores(I)=hScores(J(I));
[I J]=ismember(model.genes,arrayData.genes);
arrayScores(I)=aScores(J(I));
geneScores(I & myIsInf(geneScores))=aScores(J(I & myIsInf(geneScores)));

%Remove the genes that have no data from the model
I=ismember(model.genes,hpaData.genes) | ismember(model.genes,arrayData.genes);
model.genes(~I)=[];
model.rxnGeneMat(:,~I)=[];

%Map the genes to the HPA/array genes
[hpaExist hpaMap]=ismember(model.genes,hpaData.genes);
[arrayExist arrayMap]=ismember(model.genes,arrayData.genes);

%Set the default scores for reactions without genes
rxnScores=ones(numel(model.rxns),1)*noGeneScore;

%Loop through the reactions and calculate the scores
for i=1:numel(model.rxns)
   %Check if it has genes
   I=find(model.rxnGeneMat(i,:));
   if any(I)
       %If any of the genes exist in hpaData, then don't use arrayData
       if any(hpaExist(I))
           %At least one gene was found in HPA
           if strcmpi(multipleGeneScoring,'best')
               rxnScores(i)=max(hScores(hpaMap(I(hpaExist(I)))));
           else
               rxnScores(i)=mean(hScores(hpaMap(I(hpaExist(I)))));
           end
       else
           %Use array data
           if any(arrayExist(I))
              %At least one gene was found in the array data
               if strcmpi(multipleGeneScoring,'best')
                   rxnScores(i)=max(aScores(arrayMap(I(arrayExist(I)))));
               else
                   rxnScores(i)=mean(aScores(arrayMap(I(arrayExist(I)))));
               end
           end
       end
   end
end
end

%This is because isinf and all returns 0x1 for empty set, which gives a
%concatenation error. Do like this instead of having many if statements
function y=myIsInf(x)
	y=isinf(x);
    if isempty(y)
        y=[];
    end
end
function y=myAll(x,dim)
	y=all(x,dim);
    if isempty(y)
        y=[];
    end
end

%This is to save the mosek problem and give a hot start
function [outModel deletedRxns metProduction fValue]=runINITinServer(model,rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,noRevLoops,params,hotStart)
% runINIT
%	Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks
%
%   model           a reference model structure
%   rxnScores       a vector of scores for the reactions in the model.
%                   Positive scores are reactions to keep and negative
%                   scores are reactions to exclude (opt, default all 0.0)
%   presentMets     cell array with unique metabolite names that the model
%                   should produce (opt, default [])
%   essentialRxns   cell array of reactions that are essential and that
%                   have to be in the resulting model. This is normally
%                   used when fitting a model to task (see fitTasks) (opt,
%                   default [])
%   prodWeight      a score that determines the value of having
%                   net-production of metabolites. This is a way of having
%                   a more functional network as it provides a reason for
%                   including bad reactions for connectivity reasons. This
%                   score is for each metabolite, and the sum of these weights
%                   and the scores for the reactions is what is optimized
%                   (opt, default 0.5)
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                   This results in fewer reactions being considered
%                   dead-ends, but all reactions in the resulting model may
%                   not be able to carry flux. If this is "false" then the
%                   equality constraints are taken from model.b. If the
%                   input model lacks exchange reactions then this should
%                   probably be "true", or a large proportion of the model
%                   would be excluded for connectivity reasons
%                   (opt, default false)
%   noRevLoops      true if reversible reactions should be constrained to
%                   only carry flux in one direction. This prevents
%                   reversible reactions from being wrongly assigned as
%                   connected (the forward and backward reactions can form a
%                   loop and therefore appear connected), but it makes the 
%                   problem significantly more computationally intensive to
%                   solve (two more integer constraints per reversible reaction)
%                   (opt, default false)
%   params          parameter structure as used by getMILPParams (opt,
%                   default [])
%   hotStart        an array containing a feasible integer solution to the
%                   problem
%
%   outModel        the resulting model structure
%   deletedRxns     reactions which were deleted by the algorithm
%   metProduction   array that indicates which of the
%                   metabolites in presentMets that could be
%                   produced 
%                   -2: metabolite name not found in model
%                   -1: metabolite found, but it could not be produced
%                   1: metabolite could be produced
%   fValue          objective value (sum of (the negative of) 
%                   reaction scores for the included reactions and
%                   prodWeight*number of produced metabolites)
%
%   This function is the actual implementation of the algorithm. See
%   getINITModel for a higher-level function for model reconstruction. See 
%   PLoS Comput Biol. 2012;8(5):e1002518 for details regarding the
%   implementation.
%
%   Usage: [outModel deletedRxns metProduction fValue]=runINITinServer(model,...
%           rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,...
%           noRevLoops,params,hotStart)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    rxnScores=zeros(numel(model.rxns),1);
end
if isempty(rxnScores)
    rxnScores=zeros(numel(model.rxns),1);
end
if nargin<3
    presentMets={};
end
if isempty(presentMets)
    presentMets={};
end
presentMets=presentMets(:);
if nargin<4
    essentialRxns={};
end
if isempty(essentialRxns)
    essentialRxns={};
end
essentialRxns=essentialRxns(:);
if nargin<5
    prodWeight=0.5;
end
if isempty(prodWeight)
    prodWeight=0.5;
end
if nargin<6
    allowExcretion=false;
end
if nargin<7
    noRevLoops=false;
end
if nargin<8
    params=[];
end
if nargin<9
    hotStart = [];
end

echo=0;
if isfield(params,'printReport')
    if params.printReport==true
        echo=3;
    end
end

if numel(presentMets)~=numel(unique(presentMets))
	dispEM('Duplicate metabolite names in presentMets');
end

%Default is that the metabolites cannot be produced
if ~isempty(presentMets)
    metProduction=ones(numel(presentMets),1)*-2;
    presentMets=upper(presentMets);
    pmIndexes=find(ismember(presentMets,upper(model.metNames)));
    metProduction(pmIndexes)=-1; %Then set that they are at least found
else
    metProduction=[];
    pmIndexes=[];
end

%The model should be in the reversible format and all relevant exchange
%reactions should be open
if isfield(model,'unconstrained')
    dispEM('Exchange metabolites are still present in the model. Use simplifyModel if this is not intended',false); 
end

%The irreversible reactions that are essential must have a flux and are therefore not
%optimized for using MILP, which reduces the problem size. However, reversible
%reactions must have a flux in one direction, so they have to stay in
%the problem. The essentiality constraint on reversible reactions is
%implemented in the same manner as for reversible reactions when
%noRevLoops==true, but with the additional constraint that C ub=-1. This
%forces one of the directions to be active.
revRxns=find(model.rev~=0);
essentialReversible=find(ismember(model.rxns(revRxns),essentialRxns));
essentialRxns=intersect(essentialRxns,model.rxns(model.rev==0));

%Convert the model to irreversible
irrevModel=convertToIrrev(model);
rxnScores=[rxnScores;rxnScores(model.rev==1)];
%These are used if noRevLoops is true
if noRevLoops==true
    forwardIndexes=find(model.rev~=0);
    backwardIndexes=(numel(model.rxns)+1:numel(irrevModel.rxns))';
else
    %Then they should only be used for essential reversible reactions
    forwardIndexes=revRxns(essentialReversible);
    backwardIndexes=essentialReversible+numel(model.rxns);
end

%Get the indexes of the essential reactions and remove them from the
%scoring vector
essentialIndex=find(ismember(irrevModel.rxns,essentialRxns));
rxnScores(essentialIndex)=[];

%Go through each of the presentMets (if they exist) and modify the S matrix
%so that each reaction which produces any of them also produces a
%corresponding fake metabolite and the opposite in the reverse direction.

%This is to deal with the fact that there is no compartment info regarding
%the presentMets. This modifies the irrevModel structure, but that is fine
%since it's the model structure that is returned.
if any(pmIndexes)
    irrevModel.metNames=upper(irrevModel.metNames);
    metsToAdd.mets=strcat({'FAKEFORPM'},num2str(pmIndexes));
    metsToAdd.metNames=metsToAdd.mets;
    metsToAdd.compartments=irrevModel.comps{1};
    
    %There is no constraints on the metabolites yet, since maybe not all of
    %them could be produced
    irrevModel=addMets(irrevModel,metsToAdd);
end

%Modify the matrix
for i=1:numel(pmIndexes)
    %Get the matching mets
    I=ismember(irrevModel.metNames,presentMets(pmIndexes(i)));
    
    %Find the reactions where any of them are used.
    [crap K L]=find(irrevModel.S(I,:));
    
    %This ugly loop is to avoid problems if a metabolite occurs several
    %times in one reaction
    KK=unique(K);
    LL=zeros(numel(KK),1);
    for j=1:numel(KK)
       LL(j)=sum(L(K==KK(j)));
    end
    irrevModel.S(numel(irrevModel.mets)-numel(pmIndexes)+i,KK)=LL;
end

%Some nice to have numbers
nMets=numel(irrevModel.mets);
nRxns=numel(irrevModel.rxns);
nEssential=numel(essentialIndex);
nNonEssential=nRxns-nEssential;
nonEssentialIndex=setdiff(1:nRxns,essentialIndex);
S=irrevModel.S;

%Add so that each non-essential reaction produces one unit of a fake metabolite
temp=sparse(1:nRxns,1:nRxns,1);
temp(essentialIndex,:)=[];
S=[S;temp];

%Add another set of reactions (will be binary) which also produce these
%fake metabolites, but with a stoichiometry of 1000
temp=sparse(1:nNonEssential,1:nNonEssential,1000);
temp=[sparse(nMets,nNonEssential);temp];
S=[S temp];

%Add reactions for net-production of (real) metabolites
if prodWeight~=0
    temp=[speye(nMets-numel(pmIndexes))*-1;sparse(nNonEssential+numel(pmIndexes),nMets-numel(pmIndexes))];
    S=[S temp];
    %To keep the number of reactions added like this
    nNetProd=nMets-numel(pmIndexes);
else
    nNetProd=0;
end

%Add constraints so that reversible reactions can only be used in one
%direction. This is done by adding the fake metabolites A, B, C for each
%reversible reaction in the following manner
% forward: A + .. => ...
% backwards: B + ... => ...
% int1: C => 1000 A
% int2: C => 1000 B
% A ub=999.9
% B ub=999.9
% C lb=-1
% int1 and int2 are binary
if any(forwardIndexes)
    nRevBounds=numel(forwardIndexes);
    
    %Add the A metabolites for the forward reactions and the B
    %metabolites for the reverse reactions
    I=speye(numel(irrevModel.rxns))*-1;
    temp=[I(forwardIndexes,:);I(backwardIndexes,:)];
   
    %Padding
    temp=[temp sparse(size(temp,1),size(S,2)-numel(irrevModel.rxns))];
    
    %Add the int1 & int2 reactions that produce A and B
    temp=[temp speye(nRevBounds*2)*1000];
    
    %And add that they also consume C
    temp=[temp;[sparse(nRevBounds,size(S,2)) speye(nRevBounds)*-1 speye(nRevBounds)*-1]];
    
    %Add the new reactions and metabolites
    S=[S sparse(size(S,1),nRevBounds*2)];
    S=[S;temp];
else
    nRevBounds=0;
end

%Add so that the essential reactions must have a small flux and that the
%binary ones (and net-production reactions) may have zero flux. The
%integer reactions for reversible reactions have [0 1]
prob.blx=[irrevModel.lb;zeros(nNonEssential+nNetProd+nRevBounds*2,1)];
prob.blx(essentialIndex)=max(0.1,prob.blx(essentialIndex));

%Add so that the binary ones and net-production reactions can have at the most flux 1.0
prob.bux=[irrevModel.ub;ones(nNonEssential+nNetProd+nRevBounds*2,1)];

%Add that the fake metabolites must be produced in a small amount and that
%the A and B metabolites for reversible reactions can be [0 999.9] and C
%metabolites [-1 0]
prob.blc=[irrevModel.b(:,1);ones(nNonEssential,1);zeros(nRevBounds*2,1);ones(nRevBounds,1)*-1];

%Add that normal metabolites can be freely excreted if allowExcretion==true,
%and that the fake ones can be excreted 1000 units at most. C metabolites
%for essential reversible reactions should have an upper bound of -1.
%If noRevLoops is false, then add this constraint for all the reactions instead.
if noRevLoops==true
    revUB=zeros(nRevBounds,1);
    revUB(essentialReversible)=-1;
else
    revUB=ones(nRevBounds,1)*-1;
end
if allowExcretion==true
    metUB=inf(nMets,1);
else
    metUB=irrevModel.b(:,min(size(irrevModel.b,2),2));
end
prob.buc=[metUB;ones(nNonEssential,1)*1000;ones(nRevBounds*2,1)*999.9;revUB];

%Add objective coefficients for the binary reactions. The negative is used
%since we're minimizing. The negative is taken for the prodWeight
%as well, in order to be consistent with the syntax that positive scores
%are good
prob.c=[zeros(nRxns,1);rxnScores;ones(nNetProd,1)*prodWeight*-1;zeros(nRevBounds*2,1)];
prob.a=S;

%We still don't know which of the presentMets that can be produced. Go
%through them, force production, and see if the problem can be solved
params.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_FREE_SIMPLEX';
for i=1:numel(pmIndexes)
    prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=1;
    [crap,res] = mosekopt('minimize echo(0)', prob,getMILPParams(params));
    isFeasible=checkSolution(res);
    if ~isFeasible
        %Reset the constraint again
        prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=0;
    else
        %Metabolite produced
        metProduction(pmIndexes(i))=1;
    end
end

%Add that the binary reactions may only take integer values.
allInt=[(nRxns+1):(nRxns+nNonEssential) size(S,2)-nRevBounds*2+1:size(S,2)];
prob.ints.sub=allInt;

%Get params and provide a hot start if available
mosekParam = getMILPParams(params);

if ~isempty(hotStart)
    if size(prob.a,2) ~= max(size(hotStart))
        fprintf('WARNING: hotStart must have the same dimensions as the MILP formulation. The hot start will be ignored\n')
    else
        prob.sol.int.xx = hotStart;
    end
end

[crap,res] = mosekopt(['minimize echo(' num2str(echo) ')'], prob,mosekParam);

%I don't think that this problem can be infeasible, so this is mainly a way
%of checking the licence stuff
if ~checkSolution(res)
    dispEM('The problem is infeasible');
end

fValue=res.sol.int.pobjval;

%Get all reactions used in the irreversible model
usedRxns=(nonEssentialIndex(res.sol.int.xx(nRxns+1:nRxns+nNonEssential)<0.1))';

%Map to reversible model IDs
usedRxns=[usedRxns(usedRxns<=numel(model.rxns));revRxns(usedRxns(usedRxns>numel(model.rxns))-numel(model.rxns))];

%Then get the ones that are not used in either direction or is essential
I=true(numel(model.rxns),1);
I(usedRxns)=false;
I(essentialIndex)=false;
deletedRxns=model.rxns(I);
outModel=removeRxns(model,I,true,true);
end

%This is to parse a tab delimited file of scores in an array structure
function arrayData = parseArray(filename,tissues,celltypes,types,reliabilities)

% function parseArray
%
% The function accepts a csv delimited 2 column text file (with single line header) for which genes
% are listed in the first column and categorical levels in the second. It
% returns an arrayData structure suitable for the scoreModel RAVEN
% function.
%

fid = fopen(filename,'r');
s = textscan(fid,'%s%s','Delimiter',',','HeaderLines',1);
fclose(fid);

gene2Level_c = s{2};
nGenes = length(gene2Level_c);
arrayData.levels = unique(gene2Level_c);
nLevels = numel(arrayData.levels);
arrayData.genes = s{1};
gene2Level = zeros(nGenes,1);
for i = 1:nLevels
    gene2Level(ismember(gene2Level_c,arrayData.levels(i))) = i;
end
arrayData.gene2Level = sparse(gene2Level);

if nargin > 1
    arrayData.tissues = {tissues};
    arrayData.celltypes = {celltypes};
    arrayData.types = {types};
    arrayData.reliabilities = {reliabilities};
    arrayData.gene2Type = sparse(ones(nGenes,1));
    arrayData.gene2Reliability = sparse(ones(nGenes,1));
end
end
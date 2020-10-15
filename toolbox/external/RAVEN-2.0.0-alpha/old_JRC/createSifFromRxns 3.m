function createSifFromRxns(model,rxns,out,compartments,use_rev,fid)

% function createSifFromRxns
%
% This function creates a Cytoscape-readable .sif file for
% metabolite-gene/reaction networks. Only a GEM must be specified. If run
% on default, the resulting file is a non-compartmentalized metabolite-gene 
% non-directional network of all genes found in the model.
% The output is saved in the current directory with the name
% MetsToGenes.sif'.
%
% INPUT
% model:    RAVEN-format GEM
% rxns:     (optional) A cell array of strings containing the list of
%           reactions for which the network should be reconstructed. Default is all
%           reactions in the model. The cell array of strings must match the .rxns
%           field in the model to be included in the file.
% out:      (optional) A string, either 'rm','gm' (default) or 'rg' whether a
%           metabolite-reaction network, a metabolite-gene network, or a gene-reaction
%           network should be created.
% compartments: (optional) a logical, indicating whether metabolites in the
%               file should be considered unique if belonging to different compartments.
%               Default is false.
% use_rev:  (optional) A logical, indicating whether reversible reactions are 
%           written in the file in both directions. Default is false.
% fid:      (optional) A file ID for the output file. Default is
%           "MetsToGenes.sif"

if nargin < 6 || isempty(fid)
    fid = fopen('MetsToGenes.sif','wt');
    fprintf('\nFile saved in "MetsToGenes.sif".\n')
end
if nargin < 5 || isempty(use_rev)
    use_rev = false;
end
if nargin < 4 || isempty(compartments)
    compartments = false;
end
if nargin < 3 || isempty(out)
    out = 'gm';
end
if nargin < 2 || isempty(rxns)
    rxns = model.rxns;
end

if ~strcmp(out,{'rm';'gm';'rg'})
    fprintf('\nError: invalid out argument. It must be either "rm", "gm", or "rg".\n')
    return
end

genes     = model.genes;
mets      = model.metNames;
rev       = model.rev;
nMets     = length(mets);

if compartments
    comps       = model.comps(model.metComps); %model.compNames(str2double(model.metComps));
    metNamesC   = cellfun(@(a,b,c,d) [a,b,c,d],mets,repmat({'['},nMets,1),comps,repmat({']'},nMets,1),'uni',false);
    metNames    = metNamesC;
else
    metNames    = mets;
end

nRxns = length(rxns);
unmappedCounter = 0;
for i = 1:nRxns
    currentRxn = rxns(i);
    [mapped,indRinS] = ismember(currentRxn,model.rxns);
    if ~mapped
        unmappedCounter = unmappedCounter + 1;
        unmappedRxns(unmappedCounter) = currentRxn;
    else
        indProductsinRxn    = model.S(:,indRinS)>0;
        indSubstratesinRxn  = model.S(:,indRinS)<0;
        productsInRxn     	= metNames(indProductsinRxn);
        substratesInRxn     = metNames(indSubstratesinRxn);
        nSubs = length(substratesInRxn);
        nProd = length(productsInRxn);
        switch out
            case 'rm'
                switch use_rev                    
                    case 0
                        for j = 1:nSubs
                            fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                        end
                        for k = 1:nProd
                            fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                        end
                    case 1
                        if rev(indRinS)
                            for j = 1:nSubs
                                fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',substratesInRxn{j});
                            end
                            for k = 1:nProd
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                                fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'mr',currentRxn{:});
                            end
                        else
                            for j = 1:nSubs
                                fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                            end
                            for k = 1:nProd
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                            end
                        end
                end
            case 'gm'
                indGenesinRxn = model.rxnGeneMat(indRinS,:)~=0;
                genesInRxn    = genes(indGenesinRxn);
                nGenes        = length(genesInRxn);
                    switch use_rev                    
                        case 0
                            for l = 1:nGenes
                                for j = 1:nSubs
                                    fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',genesInRxn{l});
                                end
                                for k = 1:nProd
                                    fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'mr',genesInRxn{l});
                                end
                            end
                        case 1
                            if rev(indRinS)
                                for l = 1:nGenes
                                    for j = 1:nSubs
                                        fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',genesInRxn{l});
                                        fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'mr',substratesInRxn{j});
                                    end
                                    for k = 1:nProd
                                        fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'mr',productsInRxn{k});
                                        fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'mr',genesInRxn{l});
                                    end
                                end
                            else
                                for l = 1:nGenes
                                    for j = 1:nSubs
                                        fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',genesInRxn{l});
                                    end
                                    for k = 1:nProd
                                        fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'mr',productsInRxn{k});
                                    end
                                end
                            end
                    end
                case 'rg'
                indGenesinRxn = model.rxnGeneMat(indRinS,:)~=0;
                genesInRxn    = genes(indGenesinRxn);
                nGenes        = length(genesInRxn);
                    switch use_rev                    
                        case 0
                            for l = 1:nGenes
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',genesInRxn{l});
                            end
                        case 1
                            if rev(indRinS)
                                for l = 1:nGenes
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',genesInRxn{l});
                                    fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'mr',currentRxn{:});
                                end
                            else
                                for l = 1:nGenes
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',genesInRxn{l});
                                end
                            end
                    end
                    
        end
    end
end
fclose(fid);

if unmappedCounter > 0
    fprintf('\nThe following reactions could not be matched in the input model: ')
    disp(unmappedRxns)
    fprintf('\n')
else
    fprintf('\nAll reactions were matched to input model.\n')
end
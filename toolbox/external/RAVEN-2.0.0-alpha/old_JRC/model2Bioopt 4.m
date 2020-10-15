function SBML2Bioopt(model, BiooptFile, rxnGeneFile, sortByPathways)
%Generates a Bioopt file from a SBML model. Assumes that the model is in
%"the new" model format with respect to how gene information if stored. 
%This should be expanded on later.
%
%Primarily intended for use for the Sysbio homepage. Not to be distributed.
%
%NOTE: I assume that all characters can be used in metabolite/reaction
%names. I don't have access to Bioopt so this has to be checked.
%
%NOTE: I have to add support for objectives. Need to check the Bioopt
%manual
%
%SBML2Bioopt(SBMLfile, BiooptFile, sortByPathways)
%
%SBMLfile       A SBML file to import
%BiooptFile     The resulting Bioopt file will be saved as this
%sortByPathways True if the model should be sorted after pathways.
%               Otherwise the order in the SBML file will be used (opt, default false)
%
% Rasmus Ågren 09/08/09

if nargin<4
    sortByPathways=false;
end

[file, message]=fopen(BiooptFile,'wt');

if ~isempty(message)
    badFile=strrep(BiooptFile,'\','\\');
    fprintf(['ERROR: Could not open ' badFile ' for writing\n']);
    return;
end


if isempty(model)
	badFile=strrep(SBMLfile,'\','\\');
    fprintf(['ERROR: Could not open ' badFile ' for reading\n']);
    fclose(file);
    return;
end

%Check that no metabolites contain ':' since it's a reserved character in
%some of the tools. Print warnings and replace with '-'
badMetNames=find(~cellfun('isempty',strfind(model.metNames,':')));

if any(badMetNames)
    fprintf('WARNING: The following metabolites contains the character ":". This is a reserved character in some of the tools.\nThe metabolites have been renamed to use "-" instead\n');
    for i=1:numel(badMetNames)
       fprintf([model.metNames{badMetNames(i)} '\n']);
    end
    model.metNames=strrep(model.metNames,':','-');
end

%Get the equations for all reactions
equations=constructEquations(model);
equations=strrep(equations,'=> ','-> ');

%Print warnings for too long equations
tooLong=find(cellfun('length',equations)>=1000);

if ~isempty(tooLong)
    fprintf('WARNING: The following equations are more than 1000 characters long.\nThis is not compatible with Bioopt. All other tools function as normally.\nConsider splitting the equations into shorters fragments by introducing\nintermediate metabolites.\n\n');
    for i=1:numel(tooLong)
        fprintf([model.rxns{tooLong(i)} '\n\t' equations{tooLong(i)} '\n']);
    end
end
%Print reactions
fprintf(file,'-REACTIONS\n\n');
pad='';
for i=1:numel(model.rxns)
    fprintf(file,[pad model.rxns{i} '\t: ' equations{i} '\n']);
end

%Print constraints
fprintf(file,'\n\n-CONSTRAINTS\n\n');

%Don't print constraints for irreversible reactions if they are [0 1000] or
%reversible if they are [-1000 1000]
toPrint=find((model.rev==0 & (model.lb~=0 | model.ub~=1000)) | (model.rev==1 & (model.lb~=-1000 | model.ub~=1000)));

for i=1:numel(toPrint)
    fprintf(file,[pad model.rxns{toPrint(i)} ' [' num2str(model.lb(toPrint(i))) ', ' num2str(model.ub(toPrint(i))) ']\n']);
end

%Print unbalanced metabolites
fprintf(file,'\n\n-EXTERNAL METABOLITES\n\n');
if isfield(model,'unconstrained')
    toPrint=find(model.unconstrained);

    for i=1:numel(toPrint)
        fprintf(file,[pad model.metNames{toPrint(i)} '[' model.compNames{str2double(model.metComps{toPrint(i)})} ']\n']);
    end
end
%Prints objectives
fprintf(file,'\n\n-OBJECTIVE\n\n');
objectives=find(model.c);

if ~isempty(objectives)
    fprintf(file,'max');
    for i=1:numel(objectives)
        if i==1
            if model.c(objectives(i))>0
                fprintf(file,[' ' num2str(model.c(objectives(i))) ' ' model.rxns{objectives(i)}]);
            else
                fprintf(file,[' - ' num2str(abs(model.c(objectives(i)))) ' ' model.rxns{objectives(i)}]);
            end
        else
            if model.c(objectives(i))>0
                fprintf(file,[' + ' num2str(model.c(objectives(i))) ' ' model.rxns{objectives(i)}]);
            else
                fprintf(file,[' - ' num2str(abs(model.c(objectives(i)))) ' ' model.rxns{objectives(i)}]);
            end
        end
    end
end

%Close the file
fclose(file);

%Print the reaction-gene relationships
[file, message]=fopen(rxnGeneFile,'wt');

if ~isempty(message)
    badFile=strrep(rxnGeneFile,'\','\\');
    fprintf(['ERROR: Could not open ' badFile ' for writing\n']);
    return;
end

if isfield(model,'rxnGeneMat');
    %Get the gene encoded reactions
    geneEncoded=find(sum(model.rxnGeneMat,2));
    for i=1:numel(geneEncoded)
        fprintf(file,model.rxns{geneEncoded(i)});
        genes=find(model.rxnGeneMat(geneEncoded(i),:));
        for j=1:numel(genes)
           fprintf(file,['\t' model.genes{genes(j)}]); 
        end
        fprintf(file,'\n');
    end
else
    fprintf('ERROR: The model file contains no rxn-gene relationships. No rxn-gene interaction file will be written.\n');
    fclose(file);
    return;
end

fclose(file);
end
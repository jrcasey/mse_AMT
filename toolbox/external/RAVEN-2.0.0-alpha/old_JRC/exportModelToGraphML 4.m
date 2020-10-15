function exportModelToGraphML(model,fileName,rxnLabels,RR)
% exportModelToSIF
%   Exports a constraint-based model to a SIF file
%
%   model         a model structure
%   fileName      the filename  to export the model.
%   rxnLabels     cell array with labels for reactions.
%   metLabels     cell array with labels for metabolites (opt, default
%                 model.mets)
%
%   Usage: exportModelToSIF(model,fileName,graphType,rxnLabels,metLabels)
%
%   Rasmus Agren, 2012-12-03
%   Ibrahim 

graphType='rc';

for i=1:numel(model.mets)
     metLabels(i)={sprintf('%s[%s]',cell2mat(model.metNames(i)),model.comps{model.metComps(i)})};
end

 

if ~strcmpi(graphType,'rc') && ~strcmpi(graphType,'rr') && ~strcmpi(graphType,'cc')
   throw(MException('','The graph type is incorrect')); 
end

if numel(rxnLabels)~=numel(unique(rxnLabels))
   fprintf('WARNING: Not all reaction labels are unique\n'); 
end
if numel(metLabels)~=numel(unique(metLabels))
   fprintf('WARNING: Not all metabolite labels are unique\n'); 
end

if strcmpi(graphType,'rc')
   G=model.S;
   A=rxnLabels;
   B=metLabels;
end
if strcmpi(graphType,'rr')
   G=model.S'*model.S;
   A=rxnLabels;
   B=rxnLabels;
end
if strcmpi(graphType,'cc')
   %A metabolite is linked to all products of the reactions that it participates in
   %If G=model.S*model.S' then all connections will be double, which looks
   %weird
   irrevModel=convertToIrrev(model);
   G=sparse(numel(model.mets),numel(model.mets));
   for i=1:numel(model.mets)
      I=irrevModel.S(i,:)<0; %Get the reactions in which it is a substrate
      [J crap]=find(irrevModel.S(:,I)>0);
      G(J,i)=1;
   end
   A=metLabels;
   B=metLabels;
end

fid=fopen(fileName,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
fprintf(fid,'<graph label="Reporter Subnetwork" \n');
fprintf(fid,'xmlns:xlink="http://www.w3.org/1999/xlink" \n');
fprintf(fid,'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n');
fprintf(fid,'xmlns:cy="http://www.cytoscape.org"\n');
fprintf(fid,'xmlns="http://www.cs.rpi.edu/XGMML" \n');
fprintf(fid,'directed="1">\n');
for i=1:numel(A)
    
    I=find(ismember(RR.rxns,A(i)));
      
    fprintf(fid,'<node label="%s" id="%s">\n',cell2mat(A(i)),cell2mat(A(i)));
    if RR.express(I)==1 %Up regulated reactions
        fprintf(fid,'<graphics fill="#FF0000" width="0" outline="#000000" transparency="255"/>\n');
        fprintf(fid,' <att name="Report" type="string" value="%s"/>\n',cell2mat(RR.report(I)));
        
     elseif RR.express(I)==-1 %Down regulated reactions
        fprintf(fid,'<graphics fill="#0000A1" width="0" outline="#000000" transparency="255"/>\n');
       fprintf(fid,'<att name="Report" type="string" value="%s"/>\n',cell2mat(RR.report(I)));
     else %None regulated reactions
        fprintf(fid,'<graphics fill="#aaaaaa" width="0" outline="#000000" transparency="255"/>\n');
        fprintf(fid,'<att name="Report" type="string" alue="Not Expressed"/>\n');
     end
    fprintf(fid,'</node>\n');
end
for i=1:numel(B) %Color for Metabolites, fill
    fprintf(fid,'<node label="%s" id="%s">\n',cell2mat(B(i)),cell2mat(B(i)));
    fprintf(fid,'<graphics type="ELLIPSE" h="35.0" w="35.0" fill="#936c90" width="0" outline="#000000" transparency="255"/>\n');
    fprintf(fid,'</node>\n');
end
k=1;
 for i=1:size(G,2) %#ok<ALIGN>
     I=G(:,i)~=0;
     nodes=setdiff(B(I),A(i)); %Don't include connection to itself
     for j=1:numel(nodes)
         
       col=find(ismember(rxnLabels,A(i)));
       row=find(ismember(metLabels,nodes(j)));
       if model.rev(col)==true
            fprintf(fid,'<edge label="%d" source="%s" target="%s">\n',k,cell2mat(A(i)),cell2mat(nodes(j)));
            fprintf(fid,'<att name="weight" type="integer" value="8"/>\n');
               fprintf(fid,'<graphics width="2" fill="#1f1f1f" transparency="50"/>\n');
               fprintf(fid,'</edge>\n');
            k=k+1;
            fprintf(fid,'<edge label="%d" source="%s" target="%s">\n',k,cell2mat(nodes(j)),cell2mat(A(i)));
            fprintf(fid,'<att name="weight" type="integer" value="8"/>\n');
               fprintf(fid,'<graphics width="2" fill="#1f1f1f" transparency="50"/>\n');
               fprintf(fid,'</edge>\n');
            k=k+1;
       else
           if model.S(row,col)>0
               fprintf(fid,'<edge label="%d" source="%s" target="%s">\n',k,cell2mat(A(i)),cell2mat(nodes(j)));
               fprintf(fid,'<att name="weight" type="integer" value="8"/>\n');
               fprintf(fid,'<graphics width="2" fill="#1f1f1f" transparency="50"/>\n');
               fprintf(fid,'</edge>\n');
               k=k+1;
           else
             fprintf(fid,'<edge label="%d" source="%s" target="%s">\n',k,cell2mat(nodes(j)),cell2mat(A(i)));
             fprintf(fid,'<att name="weight" type="integer" value="8"/>\n');
             fprintf(fid,'<graphics width="2" fill="#1f1f1f" transparency="50"/>\n');
             fprintf(fid,'</edge>\n');
             k=k+1;
           end
       end
               
            
     end
     
        
%          nNodes=numel(nodes);
%          nodes(1:nNodes-1)=strcat(nodes(1:nNodes-1),{'\t'});
%          fullString=[nodes{:}];
%          fprintf(fid,[A{i} '\t' graphType '\t' fullString '\n']);
      k=k+1;
   end
fprintf(fid,'</graph>\n');
fclose(fid);


function url=drawKEGGPathwayKO(model,pathwayID)
%% temp
%model = PanGEM;
%pathwayID = 'map00010';

%%
%the aim of this function is to color the KEGG pathway based on model
% model is RAVEN model
%pathwayID is KEGG map ID such as map00010
%url link you past in web browser.
%Example
%  url=drawKEGGPathway(bif,'map00020')
%Get KO and EC from pathway ID
urlKO=sprintf('http://rest.kegg.jp/link/ko/%s',pathwayID);
% urlEC=sprintf('http://rest.kegg.jp/link/rn/%s',pathwayID);
urlRN=sprintf('http://rest.kegg.jp/link/rn/%s',pathwayID);

urlwrite(urlKO,'mapKO');
% urlwrite(urlEC,'mapEC');
urlwrite(urlRN,'mapRN');

[map mapKO]=textread('mapKO','%s%s');
[map mapRN]=textread('mapRN','%s%s');

mapKO=regexprep(mapKO,'ko:','');
mapRN=regexprep(mapRN,'rn:','');
ref=union(mapKO,mapRN);

%open EC KO map
[RN KO]=textread('kegg_ko_ec.txt','%s%s');
[KO RN]=textread('CBIOMES/DatabaseLinks/KEGG/KEGG_KO_RN.csv','%s%s','delimiter',',');

KO=regexprep(KO,'ko:','');
RN=regexprep(RN,'rn:','');


modelKO=model.genes;
k=1;
for i=1:numel(modelKO)
     ko=regexp(cell2mat(modelKO(i)),'K\d*','match');
     for j=1:numel(ko)
         newRef(k)=ko(j);
         k=k+1;
     end
     if numel(ko)==0
         ko=regexp(cell2mat(modelKO(i)),'\d*.\d*.\d*.\d*','match');
          for j=1:numel(ko)
              I =find(ismember(RN,ko));
              for kk=1:numel(I)
                  newRef(k)=RN(I(kk));
                  k=k+1;
              end
          end
     end
         
end

%find the intersection between EC and KO in the model with map EC and KO
obj=intersect(newRef,ref);

%build url
url=sprintf('http://www.kegg.jp/pathway/%s',pathwayID);
for i=1:numel(obj)
   url=sprintf('%s+%s',url,cell2mat(obj(i)));
    
end

web(url,'-browser');

function mapRN = getKEGGrxns(model,pathwayID)

%the aim of this function is to retrieve the KEGG reactions numbers for each pathway based on model
% model is RAVEN model
%pathwayID is KEGG map ID such as map00010
%url link you paste in web browser.
%Example
%  url=drawKEGGPathway(bif,'map00020')
%Get KO and EC from pathway ID
urlRN=sprintf('http://rest.kegg.jp/link/rn/%s',pathwayID);
%urlEC=sprintf('http://rest.kegg.jp/link/ec/%s',pathwayID);

urlwrite(urlRN,'mapRN');
%urlwrite(urlEC,'mapEC');

[map mapRN]=textread('mapRN','%s%s');
%[map mapEC]=textread('mapEC','%s%s');

mapRN=regexprep(mapRN,'rn:','');
%mapEC=regexprep(mapEC,'ec:','');
%ref=union(mapRN,mapEC);
% 
% %open EC RN map
% [EC KO]=textread('kegg_ko_ec.txt','%s%s');
% KO=regexprep(KO,'ko:','');
% EC=regexprep(EC,'ec:','');
% 
% 
% modelKO=model.eccodes;
% k=1;
% for i=2:numel(modelKO)
%      ko=regexp(cell2mat(modelKO(i)),'K\d*','match');
%      for j=1:numel(ko)
%          newRef(k)=ko(j);
%          k=k+1;
%      end
%      if numel(ko)==0
%          ec=regexp(cell2mat(modelKO(i)),'\d*.\d*.\d*.\d*','match');
%           for j=1:numel(ec)
%               I =find(ismember(EC,ec));
%               for kk=1:numel(I)
%                   newRef(k)=KO(I(kk));
%                   k=k+1;
%               end
%           end
%      end
%          
% end
% 
% %find the intersection between EC and KO in the model with map EC and KO
% obj=intersect(newRef,ref);
% 
% %build url
% url=sprintf('http://www.kegg.jp/pathway/%s',pathwayID);
% for i=1:numel(obj)
%    url=sprintf('%s+%s',url,cell2mat(obj(i)));
%     
% end
% 
% web(url,'-browser');

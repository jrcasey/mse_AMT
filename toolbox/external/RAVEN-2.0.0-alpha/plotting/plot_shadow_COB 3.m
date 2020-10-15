function [sneg_idx,spos_idx] = plot_shadow_COB(model,Soln,ifig,wrt)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Input : Model:RAVEN or COBRA strcuture
% % Soln  : Solution from GLPX(RAVEN or COBRA {Solutions})
% % ifig  : Figure number
% % wrt   : Logical, 'true' or 'false' : Write to a file if 'true'
% % Output:
% % sneg_idx : index of (-) shadow price for example: model.metNames(sneg_idx)
% % spos_idx : index of (+) shadow price for example: model.metNames(spos_idx)
% %
% % Written by Partho Sen 18/09/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(ifig)

subplot(1,2,1)
sneg_idx = find(Soln.y < 0);
bar(Soln.y(sneg_idx),0.5,'m'); %% negative shadow price

title({'Dual cost (- ve) Shadow price of Metabolites';'COBRA(GLPX)'},'fontsize',12);
ylabel('(Shadow price)','fontsize',12);
text(10,-1,strcat('Biomass of formation','::',num2str(Soln.f)),'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)

spos_idx = find(Soln.y > 0);

bar((Soln.y(spos_idx)),'g'); %% positive shadow price

title({'Dual cost (+ ve) Shadow price of Metabolites';'COBRA(GLPX)'},'fontsize',12);
ylabel('(Shadow price)','fontsize',12);
text(1,0.8,strcat('Biomass of formation','::',num2str(Soln.f)),'fontsize',10);
set(gca,'Xtick',[1:length(model.mets(find(Soln.y > 0)))],'Xticklabel',regexprep(model.mets(find(Soln.y > 0)),'_','-'),'fontsize',10);
rotateXLabels(gca(),30);

e=std(Soln.y(find(Soln.y > 0))).*ones(size(Soln.y(find(Soln.y > 0))));
hold on
h=errorbar(Soln.y(find(Soln.y > 0)),e);
set(h(1),'color','r')
set(h(1),'LineStyle','none');
hold off;


%%% write shadow prices to a file %%%

if(strcmp(wrt,'true'))

    disp('Write shadow price of metabolites to ./Shadow_prices_COB.xls.....')
    
             Mets = model.metNames([spos_idx;sneg_idx]);
             Y    = Soln.y([spos_idx;sneg_idx]);
            xlswrite('Shadow_prices_COB.xls',[{'Metabolites','Shadow Price_COBRA'};Mets,num2cell(Y)]);

     disp('done !')
end

end


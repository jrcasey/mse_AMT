function [sneg_idx,spos_idx] = plot_shadow(model,Soln,ifig)

figure(ifig)

% bar(Soln.x,'g');
% 
% title('Flux Distribution','fontsize',14);

subplot(1,2,1)

sneg_idx = find(Soln.y < 0);

bar(Soln.y(sneg_idx),0.5,'b'); %% negative shadow price

% heatmap((log(abs(Soln.y(sneg_idx)))));
% hold on;
% colormap(bone);
% colorbar;


% e=std(Soln.y(find(Soln.y < 0))).*ones(size(Soln.y(find(Soln.y < 0))));
% hold on
% h=errorbar(Soln.y(find(Soln.y < 0)),e);
% set(h(1),'color','r')
% set(h(1),'LineStyle','none');
% hold off;
% 
% set(gca,'Xtick',[1:length(model.mets(find(Soln.y < 0)))],'Xticklabel',regexprep(model.mets(find(Soln.y < 0)),'_','-'),'fontsize',10);
% rotateXLabels(gca(),30);
 ylabel('(Shadow price)','fontsize',12);
 title('Dual cost -- (-)Shadow price of Metabolites','fontsize',14);
% % lgd=legend(strcat('Biomass of formation','=',num2str(Soln.f)));
% % set(lgd,'color','b');
% text(15.5,-0.60,strcat('Biomass of formation','::',num2str(Soln.f)),'fontsize',10);
% grid on;

% % hold on;
% % bar(Soln.suc(find(Soln.suc)),0.4,'c') %%% corresponds to negativeshadow
% % %%price and denotes bounds of the variables or fluxes during calculating
% % %%negative shadow price


subplot(1,2,2)

spos_idx = find(Soln.y > 0);

bar((Soln.y(spos_idx)),'c'); %% positive shadow price

% heatmap(abs(log(Soln.y(spos_idx))));
% colorbar;
% colormap(cool)

% e=std(Soln.y(find(Soln.y > 0))).*ones(size(Soln.y(find(Soln.y > 0))));
% hold on
% h=errorbar(Soln.y(find(Soln.y > 0)),e);
% set(h(1),'color','k')
% set(h(1),'LineStyle','none');
% hold off;


% set(gca,'Xtick',[],'fontsize',10);
% % rotateXLabels(gca(),30);
% %legend(model.mets(find(Soln.y > 0))');
 title('Dual cost -- (+ ve) Shadow price of Metabolites','fontsize',14);
 ylabel('(Shadow price)','fontsize',12);
 text(100,1.5,strcat('Biomass of formation','::',num2str(Soln.f)),'fontsize',10);
% grid on;

% hold on;
% bar(Soln.slc(find(Soln.slc)),0.2,'g') %%% corresponds to postive shadow price

% corresponds to positive shadow
% price and denotes bounds of Dual variables corresponding to lower constraint bounds.


%%% write shadow prices to a file %%%

Mets = model.metNames(find(Soln.y));
Y = Soln.y(find(Soln.y));
%%Y=[(log(abs(Soln.y(sneg_idx))));abs(log(Soln.y(spos_idx)))];

xlswrite('Shadow_prices123_test.xls',[{'Metabolites','Shadow Price'};Mets,num2cell(Y)])

end


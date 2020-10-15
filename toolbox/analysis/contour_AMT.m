function contour_AMT(x, y, z, var, varName, varUnits, saveMe, fileName)

fig = figure('visible','off')
imagesc(x,y,z)
xlabel('Latitude')
ylabel('Depth')
colormap('jet')
hc = colorbar
ylabel(hc, strcat(varName,' [',varUnits,']'))
title(varName)
set(gca,'FontSize',20)

if saveMe
    saveas(fig,fileName,'epsc')
end

end
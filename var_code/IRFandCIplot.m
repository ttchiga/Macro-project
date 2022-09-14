
%% Plotting function that takes in the true data and %%
%%%%%% lower confidence band and the high confidence band%%%
%%%% also put the name and horizon referred to as periods  %%%


function [] = irfandCIplot(irf_estimate, lowband, highband,figname,periods)

x_plot = (1:1:periods+1);

x_fill = [x_plot, fliplr(x_plot)]; 
y_fill = [lowband',fliplr(highband')];

clf
fig = gcf;
hold on 
box on
plot(x_plot,irf_estimate,'blue','LineWidth',1.2)
fill(x_fill,y_fill,1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2);
xlim([1 periods]);
legend('estimated IRF','95% confidence band')
hold off

saveas(fig,figname);
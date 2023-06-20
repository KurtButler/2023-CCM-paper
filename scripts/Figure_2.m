%% Parameters for the plot
LINEWIDTHPARAM = 1;
FONTSIZE = 15;

%% Data generative process
% Simulate a Lorenz system
X = sample_lorenz(2000);
% Grab the x1 coordinate as an observation signal
x = xyz(X);

% Compute the autocorrelation function
[a,lags] = autocorr(x,50);

% Compute the time-delay embedding for tau = 1, 10, 40
M1 = embed(x,3,1);
M10 = embed(x,3,10);
M40 = embed(x,3,40);


%% Draw plots
figure('Position', [114 1006 1130 219])
tiledlayout(1,5,'Padding','tight','TileSpacing','compact');

nexttile;
[x,y,z] = xyz(X);
plot3(x,y,z,'k','LineWidth',LINEWIDTHPARAM)
view([1 -1 1])
xlabel('x_{1,t}','FontSize',FONTSIZE)
ylabel('x_{2,t}','FontSize',FONTSIZE)
zlabel('x_{3,t}','FontSize',FONTSIZE)
xticks([])
yticks([])
zticks([])


nexttile;
plot(lags,a,'k','LineWidth',LINEWIDTHPARAM)
hold on;
scatter([1,10,40],a(1+[1 10 40]),90,'filled')
hold off;
grid on;
xlabel('Lag','FontSize',FONTSIZE)
ylabel('ACF','FontSize',FONTSIZE)
xlim([0 50])



nexttile;
[x,y,z] = xyz(M1);
plot3(x,y,z,'k','LineWidth',LINEWIDTHPARAM)
view([1 -1 1])
xlabel('a_{t-2\tau}','FontSize',FONTSIZE)
ylabel('a_{t-\tau}','FontSize',FONTSIZE)
zlabel('a_t','FontSize',FONTSIZE)
xticks([])
yticks([])
zticks([])


nexttile;
[x,y,z] = xyz(M10);
plot3(x,y,z,'k','LineWidth',LINEWIDTHPARAM)
view([1 -1 1])
xlabel('a_{t-2\tau}','FontSize',FONTSIZE)
ylabel('a_{t-\tau}','FontSize',FONTSIZE)
zlabel('a_t','FontSize',FONTSIZE)
xticks([])
yticks([])
zticks([])


nexttile;
[x,y,z] = xyz(M40);
plot3(x,y,z,'k','LineWidth',LINEWIDTHPARAM)
view([1 -1 1])
xlabel('a_{t-2\tau}','FontSize',FONTSIZE)
ylabel('a_{t-\tau}','FontSize',FONTSIZE)
zlabel('a_t','FontSize',FONTSIZE)
xticks([])
yticks([])
zticks([])


%% Annotations
annotation("textbox",[0,0.9,0.1,0.1],'String',"(A)",'LineStyle','none','FontSize',FONTSIZE)
annotation("textbox",[0.2,0.9,0.1,0.1],'String',"(B)",'LineStyle','none','FontSize',FONTSIZE)
annotation("textbox",[0.4,0.9,0.1,0.1],'String',"(C)",'LineStyle','none','FontSize',FONTSIZE)
annotation("textbox",[0.6,0.9,0.1,0.1],'String',"(D)",'LineStyle','none','FontSize',FONTSIZE)
annotation("textbox",[0.8,0.9,0.1,0.1],'String',"(E)",'LineStyle','none','FontSize',FONTSIZE)


%% Save result
saveas(gcf,sprintf('./results/Fig2.png',date));


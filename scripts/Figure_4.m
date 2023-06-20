%% Data generative model
C = 0;
odefun = @(t,x) [-6*(x(2)+x(3))  6*(x(1)+0.2*x(2))  6*(0.2 + x(3)*(x(1)-5.7))  10*(-x(4)+x(5))  28*x(4)-x(5)-x(4)*x(6)+C*x(2)^2  x(4)*x(5)-8*x(6)/3]';
tspan = linspace(0,10,500);
X0 = [ -0.82   -0.80   -0.24    10.01    -12.19    10.70];
[t,X] = ode45(odefun,tspan,X0);

% Pick signals for CCM
X = normalize(X);
x = X(:,2);
y = X(:,1);
z = X(:,5);
t = (1:numel(t))';



%% State-Space Reconstruction
threshold = 0.5;

tauy = lag_select(y,threshold);
Qy = falsenearestneighbors(y,tauy,0.01,8);

tauz = lag_select(z,threshold);
Qz = falsenearestneighbors(z,tauz,0.01,8);

%% Convergent Cross Mapping subroutine
[CC,My,xp] = ccm(x,y,Qy,tauy);
tp  = t((numel(x)-numel(xp)) + (1:numel(xp)));
tCC = t((Qy-1)*tauy + (1:numel(CC)));


%% Plot the nice figure
[CC2,~,xp2] = ccm(x,z,Qz,tauz);
tp2  = t((numel(x)-numel(xp2)) + (1:numel(xp2)));
tCC2 = t((Qz-1)*tauz + (1:numel(CC2)));
figure('Position',[141 909 631 305])
tiledlayout(3,1,"TileSpacing","compact","Padding","tight")

nexttile
plot(t,x,t,y,t,z,'LineWidth',1);
grid on; 
grid minor;
title('','(A) Observed signals','FontSize',18)
legend('a','b','c','FontSize',15,'Location','northeast')

nexttile
plot(t,x,tp,xp,'k-.',tp2,xp2,'r--','LineWidth',1);
grid on; 
grid minor;
title('','(B) Cross predictions','FontSize',18)
legend('a','m^b\rightarrowa','m^c\rightarrowa','FontSize',15,'Location','northeast')

nexttile
plot(tCC,CC,'k',tCC2,CC2,'r','LineWidth',1);
grid on; 
grid minor;
ylim([0,1])
xlim([0,tCC(end)]);
title('','(C) Cross map skill','FontSize',18)
legend('a\Rightarrowb','a\Rightarrowc','FontSize',15,'Location','northeast')




%% Save result
saveas(gcf,sprintf('../results/Fig4.png',date));



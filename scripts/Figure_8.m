%% Code to reproduce Figure 8 of the paper
NoIter = 1000; % Number of trajectories (or realizations) per plot



ctr = 0;
for C = [0,1,3]
    %% Update counter for iteration no.
    % This is just to track which plot to write to
    ctr = ctr+1;
    
    %% Loop
    CCMmatxy = zeros(NoIter,1);  % CCM statistic for testing x-->y
    CCMmatyx = zeros(NoIter,1);  % CCM statistic for testing x<--y
    for iter = 1:NoIter
        if ~mod(iter,25); fprintf('.'); end
        %% Generate signals
        odefun = @(t,x) [-6*(x(2)+x(3))  6*(x(1)+0.2*x(2))  6*(0.2 + x(3)*(x(1)-5.7))  10*(-x(4)+x(5))  28*x(4)-x(5)-x(4)*x(6)+C*x(2)^2  x(4)*x(5)-8*x(6)/3]';
        tspan = linspace(0,50,2000);
        X0 = [0 0 0.4 0.3 0.3 0.3];
        X0 = X0 + 1e-2*randn(size(X0)); % Randomizer
        [t,X] = ode45(odefun,tspan,X0);

        % Pick signals for CCM
        X = normalize(X);
        x = X(:,2);
        y = X(:,5);

        t = t(1000:end);
        x = x(1000:end);
        y = y(1000:end);

        %% SSR
        threshold = 0.8;

        taux = lag_select(x,threshold);
        Qx = 4;

        tauy = lag_select(y,threshold);
        Qy = 8;

        %% CCM
        CCxy = ccm(x,y,Qy,tauy);
        CCyx = ccm(y,x,Qx,taux);

        % Store result
        CCMmatxy(iter) = CCxy(end);
        CCMmatyx(iter) = CCyx(end);
    end
    
    %% Add result to plot
    if ctr==1
        figure(21)
        tiledlayout(3,1,'Padding','tight','TileSpacing','compact');
    else
        figure(21)
    end
    
    nexttile 
    edges = linspace(0,1,21);
    histogram(CCMmatxy,edges,'EdgeAlpha',0.1)
    hold on;
    histogram(CCMmatyx,edges,'EdgeAlpha',0.1)
    hold off;
    grid on;
    legend('a\rightarrowb','a\leftarrowb','FontSize',13);
    title(sprintf('CCM results for C=%d',C),'FontSize',15);
    xlim([0,1])
    xlabel('Final convergence coefficient','FontSize',13)
end


%% Save result
saveas(gcf,sprintf('./results/Fig8.png',date));



%% Data generative model
C = 0;
odefun = @(t,x) [-6*(x(2)+x(3))  6*(x(1)+0.2*x(2))  6*(0.2 + x(3)*(x(1)-5.7))  10*(-x(4)+x(5))  28*x(4)-x(5)-x(4)*x(6)+C*x(2)^2  x(4)*x(5)-8*x(6)/3]';
tspan = linspace(0,15,500);
X0 = [ -0.82   -0.80   -0.24    10.01    -12.19    10.70];
[t,X] = ode45(odefun,tspan,X0);


% Observe signals
% X = normalize(X);
a = X(:,4);
N = size(a,1);



%% SSR parameter selection
threshold = 0.5;
taua = lag_select(a,threshold);
Qa = falsenearestneighbors(a,taua,0.01,8);

%% RF Analysis
% This parameter controls how much correlation we require to define a "neighbor" in the kernel-based recurrence test
corr_cutoff = 0.5;

% Initialize some guys
tRa = (Qa-1)*taua+1:N;
VectRFa = 0*tRa;

% GP Model stuff
Mmx = embed(a,Qa+1,taua);
N0 = round(size(Mmx,1)*0.8);
Mx = Mmx((1:N0),1:Qa);
mx = Mmx((1:N0),Qa+1);
gpx = fitrgp(Mx,mx);
sl = gpx.KernelInformation.KernelParameters(1); % SE kernel length scale
sf = gpx.KernelInformation.KernelParameters(2); % SE kernel amplitude
kx = @(u,v) sf.^2*exp(-0.5*pdist2(u/sl,v/sl).^2); % Kernel function

%% Recurrence analysis
% x
Mx = embed(a,Qa,taua); % Shadow manifold
Dx = pdist2(Mx,Mx); % Distance matrix
Kx = kx(Mx,Mx); % Covariance matrix
KRx= sqrt(diag(Kx)).\(Kx./sqrt(diag(Kx))); % Correlation matrix
mask = pdist2((1:size(Kx,1))',(1:size(Kx,1))')<=2*taua; % Mask for temporal neighbors

% Matrices to be plotted
oDx = Dx;
oKx = Kx;
oKRx = KRx;
oRx  = KRx>corr_cutoff & ~mask; % The recurrence matrix is spatial neighbors that are NOT temporal neighbors

KRx(mask) = 0;
Dx(mask)  = nan;

%% Plotting
figure(5)
tiledlayout(2,2,"TileSpacing","compact","Padding","tight")

nexttile
imagesc(oDx);
colorbar;
title('','Distance matrix $$\mathbf{D}$$','FontSize',15,'Interpreter','latex')

nexttile
imagesc(oKx);
colorbar;
title('','Covariance matrix $$\mathbf{K}$$','FontSize',15,'Interpreter','latex')

nexttile
imagesc(oKRx);
colorbar;
title('','Correlation matrix $$\tilde{\mathbf{K}}$$','FontSize',15,'Interpreter','latex')

nexttile
imagesc(oRx);
colorbar;
title('','Recurrence matrix $$\mathbf{R}$$','FontSize',15,'Interpreter','latex')




%% Save result
saveas(gcf,sprintf('./results/Fig6.png',date));


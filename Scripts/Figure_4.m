close all;

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



%% SSR
threshold = 0.5;

tauy = lag_select(y,threshold);
Qy = falsenearestneighbors(y,tauy,0.01,8);


tauz = lag_select(z,threshold);
Qz = falsenearestneighbors(z,tauz,0.01,8);

%% CCM
[CC,My,xp] = ccm(x,y,Qy,tauy);
tp  = t((numel(x)-numel(xp)) + (1:numel(xp)));
tCC = t((Qy-1)*tauy + (1:numel(CC)));

%% Shadow manifold
% plotembed(y,Qy,tauy);
gcf;
ans.Position= [184 235 604 533];

%% Plot the nice figure
% figure('Position',[141 909 631 305])
% tiledlayout(3,1,"TileSpacing","compact","Padding","tight")
% 
% nexttile
% plot(t,x,t,y,'LineWidth',1);
% grid on;
% grid minor;
% title('','(A) Observed signals','FontSize',15)
% 
% nexttile
% plot(t,x,tp,xp,'k-.','LineWidth',1);
% grid on;
% grid minor;
% title('','(B) Cross predictions','FontSize',15)
% 
% nexttile
% plot(tCC,CC,'k','LineWidth',1);
% ACF = autocorr(x); hold on; plot(tCC, 0*tCC + ACF(tauy),'r'); hold off;
% grid on;
% grid minor;
% ylim([0,1])
% xlim([0,tCC(end)]);
% title('','(C) Cross map skill','FontSize',15)


%% Plot the nice figure (version 2)
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
% ACF = autocorr(x); hold on; plot(tCC, 0*tCC + ACF(tauy),'r'); hold off;
grid on; 
grid minor;
ylim([0,1])
xlim([0,tCC(end)]);
title('','(C) Cross map skill','FontSize',18)
legend('a\Rightarrowb','a\Rightarrowc','FontSize',15,'Location','northeast')


%% Functions
function tau = lag_select(x,theta)
ACF = autocorr(x);
if all(ACF>=theta)
    ACF = autocorr(x,min(numel(x)-1, 1000));
end
tau = find(ACF<theta,1) - 1;
end


function s = surrogate(x,K)
    phases = exp(2i*pi*randn(size(x,1),K));
    s = real(ifft(phases.*fft(x).*flipud(conj(phases))));
end

function [CC,My,gx] = ccm(x,y,Q,tau)
% CCM
% Computes the convergent cross-mapping (CCM) for the hypothesis that
% x is causally influencing y. 
%
% Inputs:
%   x,y are Nx1 column vectors representing signals in time
%   Q       is the embedding dimension used in state-space reconstruction (SSR)
%   tau     is the embedding delay used in SSR
%
% Outputs:
%   CC  convergence coefficients estimates in time
%   My  SSR embedding of the shadow manifold of the signal y(t) 
%   gx  sequential cross-map estimates of the signal x(t)

% Normalization of input signals
xmu = mean(x);
xsig= std(x);
x = normalize(x);
y = normalize(y);

% Step 1 - State space reconstruction (delay embedding)
My = embed(y,Q,tau);

% Step 2 - Cross-mapping
n0 = (Q-1)*tau; % Time at which we have our first SSR vector
N = size(My,1); % Look at the first few guys 
K = Q+1;        % No. of points in a putative nbhd
gx = zeros(N,size(x,2)); % cross-map estimate of the time-series x(t)
for n = 1:N
    gx(n,:) = simplex(My(1:n-1,:),x(n0:n0+n-1,:), My(n,:));
end


% Save the cross-map result in the original units and timescale,
% since I provide gx(t) as an output
gx = gx*xsig + xmu;

% Step 3 - Convergence coefficients
CC = NaN(size(gx,1),1);
if size(x,2) > 1
    for t = (Q-1)*tau+1:size(gx,1)
        V1 = real(x((Q-1)*tau+(+1:t),:));
        V2 = real(gx(+1:t,:));
        CC(t) = abs(corr(  V1(:), V2(:) ,'Type','Pearson'));
    end 
else
    for t = (Q-1)*tau+1:size(gx,1)
        CC(t) = abs(corr(  real(x((Q-1)*tau+(+1:t))), real(gx(+1:t))  ,'Type','Pearson'));
    end
end
end

function Yp = simplex(M,Y,Mp)
% Simplex projection script
%   Yp = simplex(M,Y,Mp)
%
%  Inputs:
%   M  input points for training
%   Y  output points for training
%   Mp input points for prediction
%
%  Outputs:
%   Yp output predictions

if size(M,2) ~= size(Mp,2)
	error('Training and prediction inputs must have equal dimension.')
end

K = size(M,2) + 1;
Yp = zeros(size(Mp,1),size(Y,2));

for n = 1:size(Mp)
	[dists,ids] = mink( sum( (M - Mp(n,:)).^2,2 ), K);
	w = exp(-dists)/sum(exp(-dists));
	Yp(n,:) = w'*Y(ids,:);
end

end

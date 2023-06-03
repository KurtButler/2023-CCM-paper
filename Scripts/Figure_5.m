
N = 500;
NoSurro = 5;

t = (1:N)';
X = downsample(sample_lorenz(N*5),5);
x = X(:,1);
% x = awgn(x,20,'measured');
x = normalize(x);
Y = 0*x;

for subplotno = 1:2
    if subplotno==2
        x = surrogate(x);
    end

%% autoprediction
klist = 1:1:15;
tau = 5;
Q = 3;
apx = zeros(size(klist));
apY = zeros([size(klist,1),NoSurro]);

for k = klist
    apx(k) = autopred(x,Q,tau,k);
end

for s = 1:NoSurro
    for k = klist
        y = normalize( surrogate(x) );
        Y(:,s) = y;
        
    idm = (1:(N-(Q-1)*tau-k))'+(1:Q-1)*tau;
    idt = (1+(Q-1)*tau+k:N)';
        Mx = y(idm);
        mx = y(idt);
        n0 = round(size(Mx,1)*0.8);
        
        gp = fitrgp( Mx(1:n0,:), mx(1:n0,:),'BasisFunction','linear');
        mp = predict(gp, Mx(n0+1:end,:));
        mt = mx(n0+1:end,:);
        apY(k,s) = mean( (mt-mp).^2 );
        
        
        
%         apY(k,s) = autopred(y,Q,tau,k);
    end
end

%% Plot results
figure
subplot(4,1,1)
plot(x,'LineWidth',1)
grid on;
grid minor;
title('Observed signal')
xlabel('Time')

subplot(4,1,2)
plot(Y(:,1:4),'LineWidth',1)
grid on;
grid minor;
title('Surrogates')
xlabel('Time')

subplot(4,1,3)
K = max(klist);
lags = 0:K;
ax = autocorr(x,K);
aY = 0*(ax);
for s = 1:size(Y,2)
    aY(:,s) = autocorr(Y(:,s),K);
end
plot(lags,ax,'LineWidth',2)
hold on
plot(lags,aY,'k-','LineWidth',1)
hold off;
grid on;
grid minor;
legend('Observed','Surrogates','Location','best')
title('Autocorrelation function')
xlabel('Lag')
ylabel('ACF')


subplot(4,1,4);
plot(apx,'LineWidth',2);
ylim([0,1.1])

hold on;
plot(apY,'k','LineWidth',1); 
plot([tau,tau],ylim,'r','LineWidth',1)
hold off;
grid on;
grid minor;
title('Prediction MSE')
ylabel('MSE')
xlabel('Prediction length')
legend('Observed','Surrogates','Location','best')

end


%% Functions
function y = surrogate(x)
    N = size(x,1);
    if mod(N,2) % Odd N 
        m = floor(0.5*N) + 1;
        u = rand(m-1,1);
        phi = exp( 2i*pi* [0; u; -flipud(u)]);
    else
        m = floor(0.5*N);
        u = rand(m-1,1);
        phi = exp( 2i*pi* [0; u; 0; -flipud(u)]);
    end
    Y = fft(x).*phi;
    y = real(ifft(Y));
end


function MSE = autopred(x,Q,tau,k)
%PREDBENCH GPR-AR(Q) model predictability testbench
%   x    signal
%   Q    embedding dimension
%   tau  embedding delay
%   k    forecasting parameter

x = normalize(x);

M = x( (1:(size(x,1)-(Q-1)*tau-k))' + ([tau*(0:Q-1),tau*(Q-1)+k]) );
N0 = round(0.8*size(M,1));

gp = fitrgp(M(1:N0,1:end-1), M(1:N0,end),'BasisFunction','none');

m  = M(N0:end,end);
mp = predict(gp, M(N0:end,1:end-1) );

MSE = mean(abs(m - mp).^2);

end


%% Data generative process
N = 500;
t = (1:N)';
Ta = 180;
Tb = 90*sqrt(2)/1.6;
a = cos(2*pi*t/Ta);
b = cos(2*pi*t/Tb);

a = awgn(a,40,'measured');
b = awgn(b,40,'measured');

%% SSR
threshold = 0.5;

taua = 25;
Qa = 2;

taub = 11;
Qb = 2;


%% CCM
CCab = ccm(a,b,Qb,taub);
CCba = ccm(b,a,Qa,taua);

% Time indices
tCCab = (Qb-1)*taub -1 + (1:numel(CCab))';
tCCba = (Qa-1)*taua -1 + (1:numel(CCba))';

%% RF Analysis
corr_cutoff = 0.5;


tRa = (Qa-1)*taua+1:N;
tRb = (Qb-1)*taub+1:N;
VectRFa = 0*tRa;
VectRFb = 0*tRb;

% GP Models
% x
Mmx = embed(a,Qa+1,taua);
N0 = round(size(Mmx,1)*0.8);
Mx = Mmx((1:N0),1:Qa);
mx = Mmx((1:N0),Qa+1);
gpx = fitrgp(Mx,mx);
sl = gpx.KernelInformation.KernelParameters(1); % SE kernel length scale
sf = gpx.KernelInformation.KernelParameters(2); % SE kernel amplitude
kx = @(u,v) sf.^2*exp(-0.5*pdist2(u/sl,v/sl).^2); % Kernel function

% y
Mmy = embed(b,Qb+1,taub);
N0 = round(size(Mmy,1)*0.8);
My = Mmy(1:N0,1:Qb);
my = Mmy(1:N0,Qb+1);
gpy = fitrgp(My,my);
sl = gpy.KernelInformation.KernelParameters(1); % SE kernel length scale
sf = gpy.KernelInformation.KernelParameters(2); % SE kernel amplitude
ky = @(u,v) sf.^2*exp(-0.5*pdist2(u/sl,v/sl).^2); % Kernel function

%% Recurrence analysis
% x
Mx = embed(a,Qa,taua);
Dx = pdist2(Mx,Mx);
Kx = kx(Mx,Mx);
KRx= sqrt(diag(Kx)).\(Kx./sqrt(diag(Kx)));
mask = pdist2((1:size(Kx,1))',(1:size(Kx,1))')<=2*taua;
KRx(mask) = 0;
Dx(mask)  = nan;

% y
My = embed(b,Qb,taub);
Dy = pdist2(My,My);
Ky = ky(My,My);
KRy= sqrt(diag(Ky)).\(Ky./sqrt(diag(Ky)));
mask = pdist2((1:size(Ky,1))',(1:size(Ky,1))')<=2*taub;
KRy(mask) = 0;
Dy(mask)  = nan;



% Recurrence metrics
KFx = KRx>corr_cutoff;
KFy = KRy>corr_cutoff;

for l = 1:size(KFx,1)
    VectRFa(l) = mean(any(KFx(1:l,1:l)));
end
for l = 1:size(KFy,1)
    VectRFb(l) = mean(any(KFy(1:l,1:l)));
end


%% Plotting
figure(4)
tiledlayout(3,1,"TileSpacing","compact","Padding","tight")

nexttile
plot(t,a,t,b,'LineWidth',1)
grid on;
grid minor;
title('','(A) Observed signals','FontSize',15)
xlabel('Time (samples)')
ylabel('Signal amplitude')
legend('a_t','b_t','FontSize',12)


nexttile
plot([Ta,Ta],ylim, [Tb,Tb],ylim,'LineWidth',1);
hold on;
plot(tCCab,CCab,'k',tCCba,CCba,'k--','LineWidth',1)
hold off;
grid on;
grid minor;
title('','(B) Cross map skill','FontSize',15)
xlabel('Time (samples)')
ylabel('Cross map skill')
legend('','','a\Rightarrowb','a\Leftarrowb','FontSize',12)
xlim([0 500])





nexttile
plot(tRa,VectRFa,tRb,VectRFb,'LineWidth',1)
grid on;
grid minor;
title('','(C) Recurrence fractions vs time','FontSize',15)
xlabel('Time (samples)')
ylabel('RF')
legend('RF_a','RF_b','FontSize',12)




%% Save result
saveas(gcf,sprintf('./results/Fig7.png',date));



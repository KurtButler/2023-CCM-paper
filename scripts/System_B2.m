%% Parameters
NoSurro = 100;

RFparameter = 0.5;


%% Data generative process
% Bartsev example 
N = 1000;
t = (1:N)';
x = rand(size(t));
y = zeros(N,1);
for n=1:N
    if n == 1
        y(n) = rand;
    else
        y(n) = 3.4*y(n-1)*(1-y(n-1));
    end
end
a = 2 + 0.3*t./(2*250+t) + 0.1*x;
b = y;

t = t(501:1000);
a = a(501:1000);
b = b(501:1000);
N = 500;

a = normalize(a);
b = normalize(b);

sysname = "B2";

%% SSR parameters
% To make these experiments consistent with the text, we have have
% preselected the parameters to match what is written in the paper.
%   Of course, you can change the SSR parameters to experiment.
taua = 1;
taub = 1;
Qa = 4;
Qb = 4;


%% CCM
[CC,~,ga]  = ccm(a,b,Qb,taub);
[CCr,~,gb] = ccm(b,a,Qa,taua);

tCC = t(numel(t)-numel(CC) + (1:numel(CC)));
tCCr = t(numel(t)-numel(CCr) + (1:numel(CCr)));

%% Autoprediction for x-signal
k = taua; % no. of steps to predict into the future
idmnfd = (0:Qa-1)*taua + (1:(N-(Qa-1)*taua-k))';
idtrgt = (1+(Qa-1)*taua+k:N)';

Ma = zeros(size(idmnfd));
Ma(:) = a(idmnfd);
ma    = a(idtrgt);
n0 = round(size(Ma,1)*0.8);
gp = fitrgp( Ma(1:n0,:), ma(1:n0,:),'BasisFunction','linear');
r1 = sqrt(2*gp.KernelInformation.KernelParameters(1).^2*log(1/RFparameter));
mp = predict(gp, Ma(n0+1:end,:));
mt = ma(n0+1:end,:);
MSEa = mean( (mt-mp).^2 );


%% Autoprediction for y-signal
k = taub; % no. of steps to predict into the future
idmnfd = (0:Qb-1)*taub + (1:(N-(Qb-1)*taub-k))';
idtrgt = (1+(Qb-1)*taub+k:N)';

Mb = zeros(size(idmnfd));
Mb(:) = b(idmnfd);
mb    = b(idtrgt);
n0 = round(size(Mb,1)*0.8);
gp = fitrgp( Mb(1:n0,:), mb(1:n0,:),'BasisFunction','linear');
r2 = sqrt(2*gp.KernelInformation.KernelParameters(1).^2*log(1/RFparameter));
mp = predict(gp, Mb(n0+1:end,:));
mt = mb(n0+1:end,:);
MSEb = mean( (mt-mp).^2 );

%% Surrogate predictability test
surroMSEa = zeros(NoSurro,1);
surroMSEb = zeros(NoSurro,1);
asurro = zeros(size(a,1),NoSurro);
bsurro = zeros(size(b,1),NoSurro);
for s = 1:NoSurro
    asurro(:,s) = surrogate(a);
    bsurro(:,s) = surrogate(b);
end

for s = 1:NoSurro
    k = taua;
    idmnfd = (0:Qa-1)*taua + (1:(N-(Qa-1)*taua-k))';
    idtrgt = (1+(Qa-1)*taua+k:N)';
    
    Ma = zeros(size(idmnfd));
    Ma(:) = asurro(idmnfd,s);
    ma = asurro(idtrgt,s);
    n0 = round(size(Ma,1)*0.8);
    gp = fitrgp( Ma(1:n0,:), ma(1:n0,:));
    mp = predict(gp, Ma(n0+1:end,:));
    mt = ma(n0+1:end,:);
    surroMSEa(s) = mean( (mt-mp).^2 );
    
    k = taub;
    idmnfd = (0:Qb-1)*taub + (1:(N-(Qb-1)*taub-k))';
    idtrgt = (1+(Qb-1)*taub+k:N)';

    Ma = zeros(size(idmnfd));
    Ma(:) = bsurro(idmnfd,s);
    ma = bsurro(idtrgt,s);
    n0 = round(size(Ma,1)*0.8);
    gp = fitrgp( Ma(1:n0,:), ma(1:n0,:));
    mp = predict(gp, Ma(n0+1:end,:));
    mt = ma(n0+1:end,:);
    surroMSEb(s) = mean( (mt-mp).^2 );
end


%% Recurrence
Ma = embed(a,Qa,taua);
Mb = embed(b,Qb,taub);
ra = rectimes(Ma,r1,max(taua,Qa));
rb = rectimes(Mb,r2,max(taub,Qb));


%% Text to Command Window
AFa = mean(surroMSEa > MSEa);
AFb = mean(surroMSEb > MSEb);
RFa = ra(end);
RFb = rb(end);

switch (CC(end)>0.5) + 2*(CCr(end)>0.5)
    case 0
        ccmstring = "x_|_y";
    case 1
        ccmstring = "x =>y";
    case 2
        ccmstring = "x<= y";
    case 3
        ccmstring = "x<=>y";
end


fprintf('%s\t%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n',sysname,ccmstring,AFa,RFa,AFb,RFb)


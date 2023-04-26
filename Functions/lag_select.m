function tau = lag_select(x,theta)
% Select tau for SSR using the autocorrelation function
ACF = autocorr(x);
if all(ACF>=theta)
    ACF = autocorr(x,min(numel(x)-1, 100));
end
tau = find(ACF<theta,1) - 1;
end
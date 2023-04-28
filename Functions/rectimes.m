%% Functions
function rx = rectimes(Mx,r,tau)
Dx = pdist2(Mx,Mx);
tx = (1:size(Mx,1))';
Tx = pdist2(tx,tx);

% [FD,XD] = ecdf(pdist(Mx));
% r = XD(find(FD>1/sqrt(size(Mx,1)),1));
% r = norm(mean(Mx(1:end-tau,:) - Mx(1+tau:end,:)));

Rx = Dx<r & Tx> tau;

% figure, imagesc(Rx)

rx = zeros(size(Mx(:,1)));
for t = 1:numel(rx)
    rx(t) = mean(any(Rx(1:t,1:t),2));
end
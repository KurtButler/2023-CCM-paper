function Rx = recmat(Mx,r,tau)
Dx = pdist2(Mx,Mx);
tx = (1:size(Mx,1))';
Tx = pdist2(tx,tx);

% [FD,XD] = ecdf(pdist(Mx));
% r = XD(find(FD>1/sqrt(size(Mx,1)),1));
% r = norm(mean(Mx(1:end-tau,:) - Mx(1+tau:end,:)));

Rx = Dx<r & Tx> tau;
end
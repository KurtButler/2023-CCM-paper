function rx = rectimes(Mx,r,tau)
Dx = pdist2(Mx,Mx);
tx = (1:size(Mx,1))';
Tx = pdist2(tx,tx);


Rx = Dx<r & Tx> tau;


rx = zeros(size(Mx(:,1)));
for t = 1:numel(rx)
    rx(t) = mean(any(Rx(1:t,1:t),2));
end

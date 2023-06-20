function Rx = recmat(Mx,r,tau)
% This script produces the recurrence matrix with entries in {0,1}
%
% Inputs:
%   Mx = shadow manifold, generated by embed(x,Q,tau)
%   r = distance threshold for defining neighbors in space
%   tau=temporal distance for defining neighbors in time
%
% Outputs:
%   Rx = recurrence matrix

Dx = pdist2(Mx,Mx);
tx = (1:size(Mx,1))';
Tx = pdist2(tx,tx);
Rx = Dx<r & Tx> tau;
end
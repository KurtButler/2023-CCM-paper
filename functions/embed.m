function [M,t] = embed(y,Q,tau,k)
% Wrapper script to generate time-delay embedding vectors for SSR
%
% Inputs:
%   y = your signal as a column vector
%   Q = embedding dimension
%   tau=embedding lag
%   k = forecast length, if you're doing that
%
% Outputs:
%   M = the shadow manifold, as an L*Q matrix where L=1+(Q-1)*tau
%   t = forecast targets, if you're doing that

if nargin < 4
    k = 0;
end
N = size(y,1);
idm = (1:tau:1+(Q-1)*tau) + (1:N-(Q-1)*tau-1-k)';
idt = idm(:,end) + k;
M = y(idm);
t = y(idt);
end
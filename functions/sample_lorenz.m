function [X] = sample_lorenz(N,param,x0)
%SAMPLE_LORENZ 
% Samples N points along a trajectory in a Lorenz system with parameters
% given by param. Points are sampled using a Runge-Kutta method.
%   param is a vector such that param(1:3) correspond to parameters of
%   the Lorenz system, and param(4) is the sampling period of the
%   simulation method. 
%   SYNTAX          sample_lorenz(N,param,x0)
%
% Inputs:
%   N = number of points in the output signal
%   param = 4x1 vector containing [a,b,c,dt]
%   x0 = initial point for the Lorenz system
%
% Outputs:
%   X = Nx3 matrix of points along the simulated trajectory

if nargin < 2
    param = [10, 8/3, 28,0.01];
    x0 = [5 0 10];
elseif nargin < 3
    x0 = [5 0 10];
end
a = param(1);
b = param(2);
c = param(3);
dt = param(4);
X = zeros(N,3);
X(1,:) = x0(:);%+ randn(1,3); % initial state    
f = @(x) [a*(x(2)-x(1)),   x(1)*(c-x(3))-x(2),    x(1)*x(2)-b*x(3)]; % velocity function
for n=1:N-1
    x = X(n,:);
    k1 = f(x);
    k2 = f(  x+dt*k1/2);
    k3 = f(  x+dt*k2/2);
    k4 = f(  x+dt*k3);
    
    X(n+1,:) = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    if norm(x) > 100
        break;
    end
end
end


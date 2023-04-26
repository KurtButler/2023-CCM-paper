function [Q] = falsenearestneighbors(y,tau,FNNtol,Qmax)
%FALSENEARESTNEIGHBORS Select embedding dimension using false nearest
%neighbors criterion
% falsenearestneighbors(y,tau,FNNtol,Qmax)
if nargin < 3
    Qmax = 10;
end

rho = 17;
Q = 1;
FNNflag = false;
while ~FNNflag
    Q = Q + 1;
    
    if Q > Qmax
        fprintf('FNN algorithm failed to converge. FNN=%0.2f\n Forcing Q=%d.\n',mean(FNN),Qmax)
        Q = Qmax;
        break;
    end
    
    M1 = embed(y,Q,tau);
    M2 = embed(y,Q+1,tau);
    % Make sure that these guys are the same size
    M1 = M1(1:size(M2,1),:);
    FNN = zeros(size(M1,1),1);
    for n = 1:size(M1,1)
        [~,id] = mink( vecnorm(M1-M1(n,:),2,2) ,2);
        Rd = norm(M1(id(2),:)-M1(n,:),2)/sqrt(Q);
        FNN(n) = norm(M2(n,:) - M2(id(2),:),2) > rho*Rd;
    end
    
    if mean(FNN) < FNNtol
        FNNflag = true;
    end
end
end


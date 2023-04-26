function [x,y,z] = xyz(X)
%XYZ returns the first 3 columns of a matrix as a way to retrive separate
%x,y,z coordinates quickly.
x = X(:,1);
y = X(:,2);
if size(X,2) > 2
    z = X(:,3);
else
    z = zeros(size(x));
end
end


function [a,b,dx,dy] = similarityTransform(imageData,objectData)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p = size(imageData,1);
n = 2*p;
u = 4;

A(n,u) = zeros;
l(n,1) = zeros;

for i = 1:p
    A(2*i-1,1) = imageData(i,2);
    A(2*i-1,2) = -1*imageData(i,3);

    A(2*i,1) = imageData(i,3);
    A(2*i,2) = imageData(i,2);

    A(2*i-1,3) = 1;
    A(2*i,4) = 1;

    l(2*i-1) = objectData(i,2);
    l(2*i) = objectData(i,3);
end

x = inv(A'*A)*(A')*l;

a = x(1,1);
b = x(2,1);
dx = x(3,1);
dy = x(4,1);

end
function [a,b,dx,dy] = similarityTransform(imageData,objectData)
    % Summary:  performs a similarity transform from image to object plane
    %
    % Input
    %       imageData:  the image space coordinates of the points
    %       objectData: the object space coordinates of the points
    %       rmse:       the RMS error of the image space coordinates
    % Output
    %       a,b,dx,dy:       the similarity transform solved parameters
n = 2*p;
u = 4;

A(n,u) = zeros;
l(n,1) = zeros;

%Creates a and observation matrices
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

%Runs linear adjustment
x = inv(A'*A)*(A')*l;

%Outputs results
a = x(1,1);
b = x(2,1);
dx = x(3,1);
dy = x(4,1);

end
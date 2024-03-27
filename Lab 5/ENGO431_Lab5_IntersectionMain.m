clc
clear
close all


% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Left Image), (3) Y-Coordinate (Left Image), 
% (4) X-Coordinate (Right Image), (5) Y-Coordinate (Right Image)
intersectionImageCoords = load("intersectionImageCoords.txt");
cVal = 152.15;
c = 153.358;
% Format: 
% (1) Xc, (2) Yc, (3) Zc, (4) w, (5) p, (6) k
intersectionEOPs = load("intersectionEOPs.txt");
intermediate1 = intersectionEOPs(1:3,:);
intermediate2 = intersectionEOPs(4:6,:) * pi/180;
EOP = [intermediate1; intermediate2];

dataleft = intersectionImageCoords(:,2:3);
dataright = intersectionImageCoords(:,4:5);

xhat(size(ImageCoords,1),3) = zeros;




for i = size(intersectionImageCoords,1)
    [xhat, residuals,Rx,M,t,scale,RValues] = performCollinearityLSA(cVal,dataright(i,:), EOP,dataleft(i,:));
    xhat(i,:) = xhat;
end




clc
clear
close all


%TEST VALIDATION DATA USED FOR THIS ADJUSTMENT

% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Left Image), (3) Y-Coordinate (Left Image), 
% (4) X-Coordinate (Right Image), (5) Y-Coordinate (Right Image)
intersectionImageCoords = load("intersectionImageCoords.txt");
cVal = 152.15;

% Format: 
% (1) Xc, (2) Yc, (3) Zc, (4) w, (5) p, (6) k
intersectionEOPs = load("intersectionEOPs.txt");
intermediate1 = intersectionEOPs(1:3,:);
intermediate2 = intersectionEOPs(4:6,:) * pi/180;
EOP = [intermediate1; intermediate2];

dataleft = intersectionImageCoords(:,2:3);
dataright = intersectionImageCoords(:,4:5);

xhatVal(size(intersectionImageCoords,1),3) = zeros;



 
for i = 1:size(intersectionImageCoords,1)

    [xhat1] = performCollinearityLSA(cVal,dataright(i,:), EOP,dataleft(i,:));
    xhat1
    xhatVal(i,:) = xhat1;
end



%REAL DATA USED FOR THIS ADJUSTMENT





intersectionCorrectedCoords = load("Ourdata.txt");
c = 153.358;

TrueintersectionEOPs = load("EOPsForIntersectionMain.txt");
intermediateTrue1 = TrueintersectionEOPs(1:3,2:3);
intermediateTrue2 = TrueintersectionEOPs(4:6,2:3);
EOPs = [intermediateTrue1; intermediateTrue2];

dataleftTrue = intersectionCorrectedCoords(:,2:3);
datarightTrue = intersectionCorrectedCoords(:,4:5);

xhat(size(intersectionCorrectedCoords,1),3) = zeros;




for i = 1:size(intersectionCorrectedCoords,1)

    [xhatTrue, rValues] = performCollinearityLSA(c,datarightTrue(i,:), EOPs,dataleftTrue(i,:));
    xhat(i, :) = xhatTrue;
    rVal(i, :) = rValues';
end

rSum = sum(sum(rVal));

meanHeight = mean(xhat(:,3));
residuals = xhat(:,3) - meanHeight;
rms(residuals)


zero(size(residuals,1),1) = zeros;

quiver(xhat(:,1),xhat(:,2),zero,residuals(:,1),0.5)
title("Vertical Residuals of Tennis Court")
xlabel("X Object Coordinate (m)")
ylabel("Y Object Coordinate (m)")

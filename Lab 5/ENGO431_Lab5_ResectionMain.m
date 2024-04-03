clc
clear
close all
format long g

%--------------------------------Validation--------------------------------
% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Image Space), (3) Y-Coordinate (Image Space), 
% (4) X-Coordinate (Object Space), (5) Y-Coordinate (Object Space), (6) Z-Coordinate (Object Space)
resectionImageAndObjectCoords = load("resectionImageAndObjectCoords.txt");

validationImgCoords = resectionImageAndObjectCoords(:, 2:3);
validationObjCoords = resectionImageAndObjectCoords(:, 4:6);

validationImgCoords = [resectionImageAndObjectCoords(:, 1) validationImgCoords];
validationObjCoords = [resectionImageAndObjectCoords(:, 1) validationObjCoords];

% Format: 
% (1) Xc, (2) Yc, (3) Zc, (4) w, (5) p, (6) k
validationResectionEOPs = load("resectionEOPs.txt");


%--------------------------------Main Data---------------------------------
% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Left Image), (3) Y-Coordinate (Left Image), 
% (4) X-Coordinate (Right Image), (5) Y-Coordinate (Right Image)
imageCoordsMain = load("imageCoordsMain.txt");

leftImageCoords = [imageCoordsMain(:, 1), imageCoordsMain(:, 2:3)];
rightImageCoords = [imageCoordsMain(:, 1), imageCoordsMain(:, 4:5)];


% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Object Space), (3) Y-Coordinate (Object Space), (4) Z-Coordinate (Object Space)
objectCoordsMain = load("objectCoordsMain.txt");


%--------------------------------Constants---------------------------------
validationC = 152.15;   % in mm
mainC = 153.358;        % in mm

validationImgPointPrecision = 10E-5;    % in m
mainImgPointPrecision = 4E-6;           % in m

validationRMSE = 0.015;     % in mm
mainRMSE = 0.004;           % in mm

imgScale = 4900;

fidcornerMain = [106,106];
rmaxMain = norm(fidcornerMain);

%--------------------------------Running Adjustment---------------------------------

disp("Left Resect")
[xhatLeft, residualsLeft,RxLeft,RValuesLeft,xhatLeftSd] = performSinglePhotoResection(leftImageCoords,objectCoordsMain,mainRMSE,mainC,imgScale,rmaxMain);
disp("Right Resect")
[xhatRight, residualsRight,RxRight,RValuesRight,xhatRightSd] = performSinglePhotoResection(rightImageCoords,objectCoordsMain,mainRMSE,mainC,imgScale,rmaxMain);
disp("Val Resect")
[xhatVal, residualsVal,RxVal,RValuesVal,xhatValSd] = performSinglePhotoResection(validationImgCoords,validationObjCoords,validationRMSE,validationC,7800,rmaxMain);

toOutput = [xhatLeft,xhatRight];
writematrix(toOutput,"EOPsForIntersectionMainTESTINGONLY.txt")
clc
clear
close all

% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Image Space), (3) Y-Coordinate (Image Space), 
% (4) X-Coordinate (Object Space), (5) Y-Coordinate (Object Space), (6) Z-Coordinate (Object Space)
resectionImageAndObjectCoords = load("resectionImageAndObjectCoords.txt");

validationImgCoords = resectionImageAndObjectCoords(:, 2:3);
validationObjCoords = resectionImageAndObjectCoords(:, 4:6);

% Format: 
% (1) Xc, (2) Yc, (3) Zc, (4) w, (5) p, (6) k
validationResectionEOPs = load("resectionEOPs.txt");


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




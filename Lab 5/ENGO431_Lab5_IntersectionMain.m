clc
clear
close all


% Format: 
% (1) Point Number, 
% (2) X-Coordinate (Left Image), (3) Y-Coordinate (Left Image), 
% (4) X-Coordinate (Right Image), (5) Y-Coordinate (Right Image)
intersectionImageCoords = load("intersectionImageCoords.txt");

validationLeftCoords = intersectionImageCoords(:, 2:3);
validationRightCoords = intersectionImageCoords(:, 4:5);

% Format: 
% (1) Xc, (2) Yc, (3) Zc, (4) w, (5) p, (6) k
intersectionEOPs = load("intersectionEOPs.txt");

validationLeftEOPs = intersectionEOPs(:, 1);
validationRightEOPs = intersectionEOPs(:, 2);


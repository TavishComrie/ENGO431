clc
clear
close all

checkdata = load('absolute orientation input data.txt');
tieModelCoords = load('tieModelCoords.txt');
controlModelCoords = load('controlModelCoords.txt');
checkModelCoords = load('checkModelCoords.txt');


%FORMAT: ID X_m Y_m Z_m X_o Y_o Z_o
%UNITS: model coordinates in mm, object coordinates in m

[xhatCheck, residualsCheck, RxCheck, dataPrimeCheck] = performLeastSquaresAdjustment(checkdata, 152.15);
[xhat, residuals, Rx, M,t,scale] = performLeastSquaresAdjustment(c);


checkObjectCoords = [];
controlObjectCoords = [];
tieObjectCoords = [];


%for i = 1:size(checkModelCoords(i,1))
%    checkObjectCoords(i,:) = ModelTransformation();
%end


%for i = 1:size(controlModelCoords(i,1))
%    controlObjectCoords(i,:) = ModelTransformation();
%end


%for i = 1:size(tieModelCoords(i,1))
%    tieObjectCoords(i,:) = ModelTransformation();
%end

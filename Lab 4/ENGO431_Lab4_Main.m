clc
clear
close all

checkdata = load('absolute orientation input data.txt');
tieModelCoords = load('tieModelCoords.txt');
controlModelCoords = load('controlModelCoords.txt');
checkModelCoords = load('checkModelCoords.txt');
basevector = load("baseVector");



%FORMAT: ID X_m Y_m Z_m X_o Y_o Z_o
%UNITS: model coordinates in mm, object coordinates in m

[xhat, residuals, Rx, M,t,scale] = performLeastSquaresAdjustment(c);
[xhatCheck, residualsCheck, RxCheck, Mcheck, tcheck, Scalecheck] = performLeastSquaresAdjustment(checkdata);

checkObjectCoords = [];
controlObjectCoords = [];
tieObjectCoords = [];


for i = 1:size(checkModelCoords(i,1))
    checkObjectCoords(i,:) = ModelTransformation(Scalecheck,Mcheck,tcheck,checkModelCoords(i,:));
end


for i = 1:size(controlModelCoords(i,1))
    controlObjectCoords(i,:) = ModelTransformation(Scalecontrol,Mcontrol,tcontrol,controlModelCoords(i,:));
end


for i = 1:size(tieModelCoords(i,1))
    tieObjectCoords(i,:) = ModelTransformation(Scaletie,Mtie,ttie,tieModelCoords(i,:));
end

VectorPCLeft = t;

VectorPCRight = Scale * M * basevector + t;

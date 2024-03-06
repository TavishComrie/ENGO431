clc
clear
close all

checkdata = load('absolute orientation input data.txt');
tieModelCoords = load('tieModelCoords.txt');
controlModelCoords = load('controlModelCoords.txt');
checkModelCoords = load('checkModelCoords.txt');
basevector = load("baseVector.txt");



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

num_variables = 3; 
residuals_matrix = reshape(residuals, num_variables, []);

residuals_x = residuals_matrix(1, :);
residuals_y = residuals_matrix(2, :);
residuals_z = residuals_matrix(3, :);

disp('Residuals_x:');
disp(residuals_x);

disp('Residuals_y:');
disp(residuals_y);

disp('Residuals_z:');
disp(residuals_z);


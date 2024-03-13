clc
clear
close all

format long G

checkdata = load('absolute orientation input data.txt');
tieModelCoords = load('tieModelCoords.txt');
controlModelCoords = load('controlModelCoords.txt');
checkModelCoords = load('checkModelCoords.txt');
basevector = load("baseVector.txt");
controlReal = load("controlPointsReal.txt");
checkReal = load("checkPointsReal.txt");


Mimage = load('rotationMatrixRO.txt');

mainData = [controlReal(:, 1), controlModelCoords(:, 2:4), controlReal(:, 2:4)];
basevector = load("baseVector.txt");

%FORMAT: ID X_m Y_m Z_m X_o Y_o Z_o
%UNITS: model coordinates in mm, object coordinates in m

[xhat, residuals, Rx, M,t,scale] = performLeastSquaresAdjustment(mainData);
[xhatCheck, residualsCheck, RxCheck, Mcheck, tcheck, Scalecheck] = performLeastSquaresAdjustment(checkdata);

checkObjectCoords = [];
controlObjectCoords = [];
tieObjectCoords = [];


for i = 1:size(checkModelCoords, 1)
    checkObjectCoords(i,:) = ModelTransformation(scale,M,t,checkModelCoords(i,:));
end


for i = 1:size(controlModelCoords,1)
    controlObjectCoords(i,:) = ModelTransformation(scale,M,t,controlModelCoords(i,:));
end


for i = 1:size(tieModelCoords,1)
    tieObjectCoords(i,:) = ModelTransformation(scale,M,t,tieModelCoords(i,:));
end


VectorPCLeft = t;

Moverall = Mimage * transpose (Mcheck);

w = atan2d(-Moverall(3,2),Moverall(3,3));
phi = asind(Moverall(3,1));
kappa = atan2d(-Moverall(2,1),Moverall(1,1));
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


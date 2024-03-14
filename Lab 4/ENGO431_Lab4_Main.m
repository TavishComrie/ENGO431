clc
clear
close all

format long G

validationdata = load('absolute orientation input data.txt');
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

[xhat, residuals, Rx, M,t,scale] = performLeastSquaresAdjustment(mainData,0.1);
[xhatCheck, residualsCheck, RxCheck, Mcheck, tcheck, Scalecheck] = performLeastSquaresAdjustment(validationdata,1);


checkObjectCoords = [];
controlObjectCoords = [];
tieObjectCoords = [];
disp("test")

for i = 1:size(validationdata, 1)
    validationObjectCoords(i,:) = ModelTransformation(Scalecheck,Mcheck,tcheck,validationdata(i,:));
end

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

VectorPCRight = scale * M * basevector + t;

vectorPC = [VectorPCLeft, VectorPCRight]

num_variables = 3; 
residuals_matrix = reshape(residuals, num_variables, []);

residuals_x = residuals_matrix(1, :);
residuals_y = residuals_matrix(2, :);
residuals_z = residuals_matrix(3, :);
rmseX = rms(residuals_x)
rmseY = rms(residuals_y)
rmseZ = rms(residuals_z)

[checkDifferences, objectPointAccuracies] = performAccuracyAssesment(vectorPC, controlObjectCoords, checkReal, checkObjectCoords, basevector, 153.358)

%clear all paramters and stored values for good practice
clc
clear
close all

format long G

%load in all files required
validationdata = load('absolute orientation input data.txt');
tieModelCoords = load('tieModelCoords.txt');
controlModelCoords = load('controlModelCoords.txt');
checkModelCoords = load('checkModelCoords.txt');
baseVectorValidation = load("validationBaseVector.txt");
basevector = load("baseVector.txt");
controlReal = load("controlPointsReal.txt");
checkReal = load("checkPointsReal.txt");
MimageL = eye(3, 3);
MimageR_Validation = load("validationRotationMatrixRO.txt");
MimageR = load('rotationMatrixRO.txt');


%store main data for this lab
mainData = [controlReal(:, 1), controlModelCoords(:, 2:4), controlReal(:, 2:4)];
%load the base vector needed
basevector = load("baseVector.txt");

%FORMAT: ID X_m Y_m Z_m X_o Y_o Z_o
%UNITS: model coordinates in mm, object coordinates in m

%Call the LSA file with the paramters passed to obtain the final results
%for the main data
[xhat, residuals, Rx, M,t,scale, RedundancyNumbers] = performLeastSquaresAdjustment(mainData,0.1);
%Run the LSA for the valdiation data to obtain the same final results to
%compare scripts
[xhatCheck, residualsValidation, RxCheck, MValidation, tValidation, scaleValidation, RedundancyNumbersCheck] = performLeastSquaresAdjustment(validationdata,1);

%Initalize empty vectors
checkObjectCoords = [];
controlObjectCoords = [];
tieObjectCoords = [];

%For loop to transform the validation data coords using the final
%adjustment paramters to object space
for i = 1:size(validationdata, 1)
    validationObjectCoords(i,:) = ModelTransformation(scaleValidation,MValidation,tValidation,validationdata(i,:));
end

%For loop to transform the check point data coords using the final
%adjustment parameters to object space
for i = 1:size(checkModelCoords, 1)
    checkObjectCoords(i,:) = ModelTransformation(scale,M,t,checkModelCoords(i,:));
end

%For loop to transform the control point data coords using the final
%adjustment parameters to object space
for i = 1:size(controlModelCoords,1)
    controlObjectCoords(i,:) = ModelTransformation(scale,M,t,controlModelCoords(i,:));
end

%For loop to transform the tie point data coords using the final
%adjustment parameters to object space
for i = 1:size(tieModelCoords,1)
    tieObjectCoords(i,:) = ModelTransformation(scale,M,t,tieModelCoords(i,:));
end

% Determining the persepctive center coordinates in object space (left
% image)
vectorPCLeftValidation = tValidation;
VectorPCLeft = t;

% Determining the persepctive center coordinates in object space (right
% image)
vectorPCRightValidation = scaleValidation * MValidation * baseVectorValidation + tValidation;
VectorPCRight = scale * M * basevector + t;

% Storing the pairs of perspective cetner coordinates into a vector
vectorPCValidation = [vectorPCLeftValidation, vectorPCRightValidation];
vectorPC = [VectorPCLeft, VectorPCRight];

% Computing the object to image space rotation matrix for the validation 
Moi_L_Validation = MimageL * transpose(MValidation);
Moi_R_Validation = MimageR_Validation * transpose(MValidation);

% Computing the object to image space rotation matrix
Moi_L = MimageL * transpose(M);
Moi_R = MimageR * transpose(M);

% Extracting rotation angles for the left image (validation)
w_L_Validation = atan2d(-Moi_L_Validation(3,2), Moi_L_Validation(3,3));
phi_L_Validation = asind(Moi_L_Validation(3,1));
kappa_L_Validation = atan2d(-Moi_L_Validation(2,1),Moi_L_Validation(1,1));

% Extracting rotation angles for the right image (validation)
w_R_Validation = atan2d(-Moi_R_Validation(3,2),Moi_R_Validation(3,3));
phi_R_Validation = asind(Moi_R_Validation(3,1));
kappa_R_Validation = atan2d(-Moi_R_Validation(2,1),Moi_R_Validation(1,1));

% Extracting rotation angles for left image 
w_L = atan2d(-Moi_L(3,2),Moi_L(3,3));
phi_L = asind(Moi_L(3,1));
kappa_L = atan2d(-Moi_L(2,1),Moi_L(1,1));

% Extracting rotation angles for right image 
w_R = atan2d(-Moi_R(3,2),Moi_R(3,3));
phi_R = asind(Moi_R(3,1));
kappa_R = atan2d(-Moi_R(2,1),Moi_R(1,1));

extractedAnglesValidation = [w_L_Validation, phi_L_Validation, kappa_L_Validation;
                             w_R_Validation, phi_R_Validation, kappa_R_Validation];

extractedAngles = [w_L, phi_L, kappa_L;
                   w_R, phi_R, kappa_R];


%initialize num variables to split vector of residuals into x y z 
num_variables = 3; 
%store residuals into matrix with columns rows x, y, z and columns = to numvariables
residuals_matrix = reshape(residualsValidation, num_variables, []);

%separate residuls into coordiante vectors and solve for rmse
residuals_x = residuals_matrix(1, :);
residuals_y = residuals_matrix(2, :);
residuals_z = residuals_matrix(3, :);
rmseX = rms(residuals_x);
rmseY = rms(residuals_y);
rmseZ = rms(residuals_z);


%Find coordinate differences for both the validation and the group data points
objectBaseVectorValidation = [vectorPCValidation(1, 2) - vectorPCValidation(1, 1);
                              vectorPCValidation(2, 2) - vectorPCValidation(2, 1);
                              vectorPCValidation(3, 2) - vectorPCValidation(3, 1)];

objectBaseVector = [vectorPC(1, 2) - vectorPC(1, 1);
                    vectorPC(2, 2) - vectorPC(2, 1);
                    vectorPC(3, 2) - vectorPC(3, 1)];

%Pass in object differnces and other paramters to find all accuracy
%assesment solututions for this lab
[checkDifferencesValidation, objectPointAccuraciesValidation] = performAccuracyAssesment(vectorPCValidation, validationObjectCoords, validationdata(:, 5:7), validationObjectCoords, objectBaseVectorValidation, 152.15, 1E-5);
[checkDifferences, objectPointAccuracies] = performAccuracyAssesment(vectorPC, controlObjectCoords, checkReal(:, 2:4), checkObjectCoords, objectBaseVector, 153.358, 0.004E-3);

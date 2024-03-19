clc
clear
close all

%formal xl,yl, xr,yr units of mm
datacheck = load("relative orientation input data.txt");
datacheck(:,1)=[];

tiePoints = load("tiePoints.txt");
controlPoints = load("controlPoints.txt");
checkPoints = load("checkPoints.txt");

c = 153.358;

[xhatValidation, residualsValidation, RxValidation, dataPrimeValidation, MValidation] = performLeastSquaresAdjustment(datacheck, 152.15);
[xhatTie, residualsTie, RxTie, dataPrimeTie, MTie] = performLeastSquaresAdjustment(tiePoints, c);

baseVectorValidation = [92.000; xhatValidation(1, 1); xhatValidation(2, 1)];
baseVector = [92.000; xhatTie(1, 1); xhatTie(2, 1)];

for i = 1:size(controlPoints, 1)
    vectorControl = [controlPoints(i, 3); controlPoints(i, 4); -c];
    dataPrimeControl(i, :) = transpose(MTie) * vectorControl;
end

for i = 1:size(checkPoints, 1)
    vectorCheck = [checkPoints(i, 3); checkPoints(i, 4); -c];
    dataPrimeCheck(i, :) = transpose(MTie) * vectorCheck;
end


writematrix(xhatValidation, 'ROPCheck.txt')
writematrix(residualsValidation, 'residualsCheck.txt')
writematrix(RxValidation, 'CorrCheck.txt')
writematrix(xhatTie, 'ROP.txt')
writematrix(residualsTie, 'residuals.txt')
writematrix(RxTie, 'Corr.txt')


%Puts into dd to check
xhatValidation(3:5,1) = xhatValidation(3:5,1) * 180 / pi;
xhatTie(3:5,1) = xhatTie(3:5,1) * 180 / pi;


% Performing Space Intersection
[xhatValidationSI, yPValidation] = performSpaceIntersection(datacheck, dataPrimeValidation, xhatValidation, 152.15)
[xhatTieSI, yPTie] = performSpaceIntersection(tiePoints, dataPrimeTie, xhatTie, c)
[xhatControlSI, yPControl] = performSpaceIntersection(controlPoints, dataPrimeControl, xhatTie, c)
[xhatCheckSI, yPCheck] = performSpaceIntersection(checkPoints, dataPrimeCheck, xhatTie, c)


writematrix(baseVectorValidation, "validationBaseVector.txt")
writematrix(baseVector, "baseVector.txt")
writematrix([ones(size(xhatTieSI, 1), 1), xhatTieSI], "tieModelCoords.txt")
writematrix([ones(size(xhatControlSI, 1), 1), xhatControlSI], "controlModelCoords.txt")
writematrix([ones(size(xhatCheckSI, 1), 1), xhatCheckSI], "checkModelCoords.txt")
writematrix(MValidation, "validationRotationMatrixRO.txt")
writematrix(MTie, "rotationMatrixRO.txt")
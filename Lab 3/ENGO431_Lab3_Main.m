clc
clear
close all

%formal xl,yl, xr,yr units of mm
datacheck = load("relative orientation input data.txt");
datacheck(:,1)=[];

datamain = load("uofc relative orientation input data.txt");

[xhatCheck, residualsCheck, RxCheck, dataPrimeCheck] = performLeastSquaresAdjustment(datacheck)
[xhatMain, residualsMain, RxMain, dataPrimeMain] = performLeastSquaresAdjustment(datamain)

writematrix(xhatCheck, 'ROPCheck.txt')
writematrix(residualsCheck, 'residualsCheck.txt')
writematrix(RxCheck, 'CorrCheck.txt')
writematrix(xhatMain, 'ROP.txt')
writematrix(residualsMain, 'residuals.txt')
writematrix(RxMain, 'Corr.txt')


%Puts into dd to check
xhatCheck(3:5,1) = xhatCheck(3:5,1) * 180 / pi
xhatMain(3:5,1) = xhatMain(3:5,1) * 180 / pi


clc
clear
close all

checkdata = load('absolute orientation input data.txt');

%FORMAT: ID X_m Y_m Z_m X_o Y_o Z_o
%UNITS: model coordinates in mm, object coordinates in m

[xhatCheck, residualsCheck, RxCheck,Mcheck,tcheck,scalecheck] = performLeastSquaresAdjustment(checkdata);


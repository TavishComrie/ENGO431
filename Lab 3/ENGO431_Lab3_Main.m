clc
clear
close all

%formal xl,yl, xr,yr units of mm
datacheck = load("relative orientation input data.txt");
datacheck(:,1)=[];

datamain = load("uofc relative orientation input data.txt");

%LeastSquaresAdustment(datacheck);
%LeastSquaresAdustment(datamain);
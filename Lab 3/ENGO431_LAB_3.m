clc
close all

focal_length = 153.358;

Data =        [30	106.399	90.426	24.848	81.824;
        40	18.989	93.365	-59.653	88.138;
        72	70.964	4.907	-15.581	-0.387;
        127	-0.931	-7.284	-85.407	-8.351;
        112	9.278	-92.926	-78.81	-92.62;
        50	98.681	-62.769	8.492	-68.873];

xprime = zeros(6,1);
yprime = zeros(6,1);
zprime = zeros(6,1);

for i = 1:size(Data,1)
    mvectorXYZ = [Data(i,4), Data(i,5), focal_length];
    M = M_transformation_Matrix(Xo);
    [xprime(i),yprime(i),zprime(i)] = transformation(M,mvectorXYZ);
end

function M = M_transformation_Matrix(Xnot)
    w = Xnot(1,1);
    phi = Xnot(2,1);
    kappa = Xnot(3,1);
    M = [cosd(phi)*cosd(kappa), cosd(w)*sind(kappa)+sind(w)*sind(phi)*cosd(kappa), sind(w)*sind(kappa)-cosd(w)*sind(phi)*cosd(kappa);
        -cosd(phi)*sind(kappa), cosd(w)*cosd(kappa)-sind(w)*sind(phi)*sind(kappa), sind(w)*cosd(kappa)+cosd(w)*sind(phi)*sind(kappa);
        sind(phi), -sind(w)*cosd(phi), cosd(w)*cosd(phi)];
end


function [xprime,yprime,zprime] = transformation(M, vector_x_y_z)
    vector_intermediate = transpose(M)* transpose(vector_x_y_z);
    xprime = vector_intermediate(1,1);
    yprime = vector_intermediate(2,1);
    zprime = vector_intermediate(3,1);
end
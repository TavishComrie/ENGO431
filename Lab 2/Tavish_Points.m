clc
close all
clear 

% %%Lab 2 Affine Transformation

%Entering Fiducial Coordinates in pixels
fids = [1346.25	-19284.75	1345.75	-19286
19161	-1470.25	19174.75	-1483.25
1345	-1471.75	1358.25	-1471.75
19162.5	-19283.75	19162.25	-19297.75
842	-10378	848.75	-10378.25
19665.25	-10376.5	19672.25	-10391.25
10254	-966.25	10267	-973
10254.25	-19788.25	10253.25	-19795.5
];

%Entering given fiducial coordinates in mm from calibration certificate
obsin = [-105.997	-105.995
106.004	106.008
-106	106.009
106.012	-105.995
-112	0.007
112.006	0.007
0.005	112.007
0.002	-111.998
];

%Put the obs in a single column of x,y,x,y variables
obs = normalize(obsin(:,1),obsin(:,2));

%Put the fiducials for image 27 and image 28 into single columns of
%x,y,x,y... for each specific image
fids27 = normalize(fids(:,1),fids(:,2));
fids28 = normalize(fids(:,3),fids(:,4));

%Creates A matrices
A27 = makeA(fids27);
A28 = makeA(fids28);


%Solves for unknowns
xhat27 = makeXhat(A27,obs);
xhat28 = makeXhat(A28,obs);

%Solves for residuals
v27 = residuals(A27,xhat27,obs);
v28 = residuals(A28,xhat28,obs);


%Determines object parameters from a,b,c,d
[theta,Sx,Sy,skew] = calcParams(xhat27(1),xhat27(2),xhat27(4),xhat27(5));
[theta28,Sx28,Sy28,skew28] = calcParams(xhat28(1),xhat28(2),xhat28(4),xhat28(5));

%Adjusts and split x,y,x,y columbs back into two columns of x and y
%coordinates
finalcoords1 = adjust(xhat27,fids27);
pCoords27 = unnorm(finalcoords1);
finalcoords2 = adjust(xhat28,fids28);
pCoords28 = unnorm(finalcoords2);

%Splits the residuals into x and y columns
res27 = unnorm(v27);
res28 = unnorm(v28);

%Finds rms error 
rmsX27 = rms(res27(:,1))
rmsY27 = rms(res27(:,2))
rmsX28 = rms(res28(:,1))
rmsY28 = rms(res28(:,2))

%Plots all of the results
figure;
quiver(pCoords27(:,1),pCoords27(:,2),res27(:,1),res27(:,2),0.1);
xlabel('x (mm)')
ylabel('y (mm)')
title('Image Point Residuals from Image 27')
print(gcf, '5', '-dpng', '-r300');

figure;
quiver(pCoords28(:,1),pCoords28(:,2),res28(:,1),res28(:,2),0.1);
xlabel('x (mm)')
ylabel('y (mm)')
title('Image Point Residuals from Image 28')
print(gcf, '6', '-dpng', '-r300');

%Control and Check points from image 27
points27 = [11781	-1174
9616.5	-14501.25
17840.75	-18026.5
	
9460	-2292
14007.75	-9748.75
10059.25	-10882.5
11844.5	-17253.5
14686	-18205.5
];

%Control and Check points from image 28
points28 = [3726	-851
1953	-14415.5
10137.25	-17688.75
	
1411.5	-2079.5
6160.5	-9527
2274.5	-10786
4158.25	-17085
6984	-17948.5
];
%Tie points from image 27
tPoints27 = [14358.25	-4412.00
19026.75	-4494.75
9293.00	-9500.50
17891.25	-10881.00
11327.00	-14557.75
19133.25	-16965.00
];

%Tie points from image 28
tPoints28 = [6404.25	-4093.75
11113.00	-4029.25
1434.75	-9419.50
10047.50	-10555.75
3633.75	-14424.75
11400.00	-16602.75
];

%Split columns into a single column of x,y,x,y then adjust points using the
%Affine parameters and then reverted the column back to two columns. For
%image 27 control and check points.
normPoints27 = normalize(points27(:,1),points27(:,2));
adjPoints27 = adjust(xhat27,normPoints27);
Fpoints27 = unnorm(adjPoints27);

%Split columns into a single column of x,y,x,y then adjust points using the
%Affine parameters and then reverted the column back to two columns. For
%image 28 control and check points
normPoints28 = normalize(points28(:,1),points28(:,2));
adjPoints28 = adjust(xhat28,normPoints28);
Fpoints28 = unnorm(adjPoints28);

%Extracting Correction information for images 27 and 28 control and check
%points
[Xprime27,Yprime27, correctionsMatrix27] = ImagePointCorrections(Fpoints27(:,1),Fpoints27(:,2));
[Xprime28,Yprime28, correctionsMatrix28] = ImagePointCorrections(Fpoints28(:,1),Fpoints28(:,2));

%Completing the same process as above for tie points of image 27
normTPoints27 = normalize(tPoints27(:,1),tPoints27(:,2));
adjTPoints27 = adjust(xhat27,normTPoints27);
FTpoints27 = unnorm(adjTPoints27);
[XprimeT27,YprimeT27, correctionsMatrixT27] = ImagePointCorrections(FTpoints27(:,1),FTpoints27(:,2));

%Completing the same process as above for tie points of image 28
normTPoints28 = normalize(tPoints28(:,1),tPoints28(:,2));
adjTPoints28 = adjust(xhat28,normTPoints28);
FTpoints28 = unnorm(adjTPoints28);
[XprimeT28,YprimeT28, correctionsMatrixT28] = ImagePointCorrections(FTpoints28(:,1),FTpoints28(:,2));

tiePoints = [XprimeT27, YprimeT27, XprimeT28, YprimeT28];

writematrix(tiePoints, "tiePoints.txt")

%Transformation check with Example from lecture notes
%Following same steps as for the previous transformation of the lab data 
Comp_cords= [-113.767	-107.400
-43.717		-108.204
36.361		-109.132
106.408		-109.923
107.189		-39.874
37.137		-39.070
-42.919		-38.158
-102.968	-37.446
-112.052	42.714
-42.005		41.903
38.051		40.985
108.089		40.189
108.884		110.221
38.846		111.029
-41.208		111.961
-111.249	112.759];

Reseaux_Coords = [-110	-110
-40	-110
40	-110
110	-110
110	-40
40	-40
-40	-40
-100	-40
-110	40
-40	40
40	40
110	40
110	110
40	110
-40	110
-110	110];


obs2 = normalize((Comp_cords(:,1)),Comp_cords(:,2));
res = normalize(Reseaux_Coords(:,1),Reseaux_Coords(:,2));



A= makeA(obs2);

xhat = makeXhat(A,res);


v = residuals(A,xhat,res);


%[theta,Sx,Sy,skew] = calcParams(xhat(1),xhat(2),xhat(4),xhat(5))
finalcoords = adjust(xhat,obs2);
pCoords = unnorm(finalcoords);

resids = unnorm(v);

figure;
quiver(pCoords(:,1),pCoords(:,2),resids(:,1),resids(:,2),0.1)
xlabel('x (mm)')
ylabel('y (mm)')
title('Validation Image Point Residual Quiver Plot')
print(gcf, '7', '-dpng', '-r300');

% Data output
controlPoints = [Xprime27(1:3, :), Yprime27(1:3, :), Xprime28(1:3, :), Yprime28(1:3, :)];
checkPoints = [Xprime27(4:end, :), Yprime27(4:end, :), Xprime28(4:end, :), Yprime28(4:end, :)];

writematrix(controlPoints, "controlPoints.txt")
writematrix(checkPoints, "checkPoints.txt")

%Function to populate the A matrix using the formula described
% in the lab handout
function A = makeA(fids)
    A = zeros(size(fids,1),6);
    for i = 1:2:size(fids,1)
            A(i,1) = fids(i);
            A(i,2) = fids(i+1);
            A(i,3) = 1;
            A(i+1,4) = fids(i);
            A(i+1,5) = fids(i+1);
            A(i+1,6) = 1;
            
    end
end

%Function to find the adjusted parameteres using the formula described
% in the lab report
function xhat = makeXhat(A,obs)
   xhat = ((transpose(A)*A)^-1)*transpose(A)*obs;
end

%Function to find the residuals using the formula described
% in the lab report
function v = residuals(A,xhat,obs)
    v = (A*xhat)-obs;
end

%Function to transform the affine linear parameter variables into
%non-linear parameters using the equations described in the lab report
function [theta,Sx,Sy,skew] = calcParams(a,b,c,d)
    theta = atan2(c,a);
    Sx = sqrt(a^2+c^2);
    Sy = sqrt(b^2+d^2);
    skew = atan2((a*b+c*d),(a*d-b*d));
end

%Function to split two columns of a matrix into one alternating column,
%takes two columns as input and interates through to add them in the
%pattern Col1 then Col2 then Col1 and so on.
function normal = normalize(col1,col2)
    normal = zeros(size(col1,1)*2,1);
    j=1;
    
    for i = 1:2:size(col1,1)*2
        normal(i) = col1(j);
        normal(i+1) = col2(j);
    
        j= j+1;
    end
end

%Function to reverse the normalize function therefore to take one
%alternating column and revert it back to two columns. Takes on vector as
%input and iterates through each row alternating between assigning them to
%a column
function reverse = unnorm(res)
    reverse = zeros((size(res,1))/2,2);
    j = 1;
    for i = 1:2:size(res,1)
        reverse(j,1) = res(i);
        reverse(j,2) = res(i+1);

        j = j+1;
    end
end

%Function to adjust pixel coordinates into mm using the affine
%transformation parameters solved for previously in this lab. Takes the
%affine parameters and a vector of coordinates as input and returns a
%vector of adjusted coordinates.
function finalcoords = adjust(xhat,cCoords)
    finalcoords = zeros(size(cCoords,1),1);
    
    R1 = [xhat(1) xhat(2)
        xhat(4), xhat(5)];
    R2 = [xhat(3)
        xhat(6)];

    for i = 1:2:size(cCoords,1)
        R3 = [cCoords(i)
            cCoords(i+1)];
       setCoords = (R1*R3)+R2;
       finalcoords(i,1) = setCoords(1);
       finalcoords(i+1,1) = setCoords(2);
    end
end
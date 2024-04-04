clc % Clear command window
close all % Close all figures
clear % Clear workspace

% Define Image28PointsMatrix
Image28PointsMatrix = [500	7398	-3018
                       501	7615	-3010
                       502	7858	-3002
                       503	8125	-2995
                       504	8128	-3088
                       505	8133	-3309
                       506	8135	-3403
                       507	7869	-3411
                       508	7627	-3418
                       509	7408	-3425
                       510	7406	-3331
                       511	7401	-3112
                       512	7642	-3104
                       513	7885	-3096
                       514	7891	-3316
                       515	7648	-3324];

% Extract ID column from Image28PointsMatrix
id = Image28PointsMatrix(:,1);

% Define Image27PointsMatrix
Image27PointsMatrix = [500	15362	-3393
                       501	15579	-3393
                       502	15817	-3393
                       503	16080	-3393
                       504	16081	-3485
                       505	16081	-3700
                       506	16081	-3792
                       507	15819	-3791
                       508	15579	-3791
                       509	15363	-3791
                       510	15363	-3699
                       511	15362	-3485
                       512	15602	-3485
                       513	15841	-3485
                       514	15841	-3699
                       515	15602	-3699];

% Define Xhat_27 and Xhat_28
Xhat_27 = [0.0118999833280590;
           5.05059523457756e-07;
          -122.010584700925;
          -8.93651633359339e-07;
           0.0119013719065776;
           123.521156433127];

Xhat_28 = [0.0118997602707726;
          -8.61376898313572e-06;
          -122.180337677249;
           7.90966041976906e-06;
           0.0119007707207354;
           123.509763906233];

% Adjust coordinates using adjust function for both sets of points
Finalcoords28 = adjust(Xhat_28, Image28PointsMatrix);
Finalcoords27 = adjust(Xhat_27, Image27PointsMatrix);

% Perform image point corrections for both sets of points
[Xprime27, yPrime27, matrix_corr27] = ImagePointCorrections(Finalcoords27(:,1), Finalcoords27(:,2));
[Xprime28, yPrime28, matrix_corr28] = ImagePointCorrections(Finalcoords28(:,1), Finalcoords28(:,2));

% Combine data for both sets of points
data27 = [id,Xprime27, yPrime27];
data28 = [Xprime28, yPrime28];
data = [data27, data28];

% Write data to file
writematrix(data, "Ourdata.txt")

% Define adjust function to adjust coordinates
function finalcoords = adjust(xhat, PointsMatrix)
    finalcoords = zeros(size(PointsMatrix) -1);

    R1 = [xhat(1) xhat(2); xhat(4) xhat(5)];
    R2 = [xhat(3); xhat(6)];

    for i = 1:size(PointsMatrix, 1)
        R3 = [PointsMatrix(i, 2); PointsMatrix(i, 3)];
        setCoords = (R1 * R3) + R2;
        finalcoords(i, :) = setCoords';
    end
end

% Define function for image point corrections
function [xPrime,yPrime,correctionsMatrix] = ImagePointCorrections(x,y)
    % Principal point offset correction
    xp=-0.006;
    yp=0.006;
    xBar = x - xp; 
    yBar = y - yp;

    % Find other corrections
    [xrad,yrad] = findRadialLensCorrection(xBar,yBar);
    [xdec,ydec] = DecentringLensDistortion(xBar,yBar);
    [xatm,yatm] = findAtmosphericRefractionCorrection(xBar,yBar);

    % Apply corrections
    xPrime = xBar + xrad + xdec + xatm;
    yPrime = yBar + yrad + ydec + yatm;

    % Store corrections matrix
    correctionsMatrix = [xrad,yrad,xdec,ydec,xatm,yatm];
end

clc
close all
clear 
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
515	7648	-3324
];




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
515	15602	-3699
];


Xhat_27 = [0.0118999833280590;
5.05059523457756e-07;
-122.010584700925;
-8.93651633359339e-07;
0.0119013719065776;
123.521156433127];

Xhat_28 = [
    0.0118997602707726;
-8.61376898313572e-06;
-122.180337677249;
7.90966041976906e-06;
0.0119007707207354;
123.509763906233];


Finalcoords28 = adjust(Xhat_28, Image28PointsMatrix);
Finalcoords27 = adjust(Xhat_27, Image27PointsMatrix);


[Xprime27, yPrime27, matrix_corr27] = ImagePointCorrections(Finalcoords27(:,1), Finalcoords27(:,2));
[Xprime28, yPrime28, matrix_corr28] = ImagePointCorrections(Finalcoords28(:,1), Finalcoords28(:,2));


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


function [xPrime,yPrime,correctionsMatrix] = ImagePointCorrections(x,y)
%ImagePointCorrections
% x and y are provided in fiducial system. Units of mm. May be a vector
% xPrime and yPrime return corrected coordinates for principal point,
% radial lens, decentring and atmoshperic factors

    %From calibration sheet
    xp=-0.006;
    yp=0.006;
    
    %Principal point offset correction
    xBar = x - xp; 
    yBar = y - yp;
    

    %Find other corrections
    [xrad,yrad] = findRadialLensCorrection(xBar,yBar);
    [xdec,ydec] = DecentringLensDistortion(xBar,yBar);
    [xatm,yatm] = findAtmosphericRefractionCorrection(xBar,yBar);

    %Apply corrections
    xPrime = xBar + xrad + xdec + xatm;
    yPrime = yBar + yrad + ydec + yatm;

    correctionsMatrix = [xrad,yrad,xdec,ydec,xatm,yatm];
end



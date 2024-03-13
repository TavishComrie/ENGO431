function [differences, accuracies] = performAccuracyAssesment(objectCoordsPC, controlPoints, realCheckPoints, calcCheckPoints, residuals, baseVector, c)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    imagePointPrecision = 0.004E-3;
    H = mean(objectCoordsPC(:, 3));
    h = mean(controlPoints(:, 3));

    S = (H - h) / c;

    baseMagnitude = sqrt(baseVector(1, 1)^2 + baseVector(2, 1)^2 + baseVector(3, 1)^2);
    BH = baseMagnitude / (H - h);



    differencesX = realCheckPoints(:, 5:7) - calcCheckPoints(:, :);

end
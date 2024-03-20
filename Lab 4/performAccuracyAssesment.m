function [differences, accuracies] = performAccuracyAssesment(objectCoordsPC, controlPoints, realCheckPoints, calcCheckPoints, baseVector, c, imagePointPrecision)
    differences = realCheckPoints - calcCheckPoints;

    H = mean(objectCoordsPC(3, :));
    h = mean(controlPoints(:, 3));

    S = (H - h) / (c * 10^-3);

    baseMagnitude = sqrt(baseVector(1, 1)^2 + baseVector(2, 1)^2 + baseVector(3, 1)^2);
    BH = baseMagnitude / (H - h)

    accuracyX = S * imagePointPrecision;
    accuracyY = accuracyX;
    accuracyZ = ((sqrt(2) * S) / BH) * imagePointPrecision;

    accuracies = [accuracyX; accuracyY; accuracyZ];
end
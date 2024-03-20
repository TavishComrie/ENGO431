function [differences, accuracies] = performAccuracyAssesment(objectCoordsPC, controlPoints, realCheckPoints, calcCheckPoints, baseVector, c, imagePointPrecision)
    % Summary: the performAccuracyAssesment function returns the 
    % 
    % Input
    %       objectCoordsPC:         the object space coordinates of the perspective
    %                               centers 
    %       controlPoints:          the given control points in object space
    %       realCheckPoints:        the given check points in object space
    %       calcCheckPoints:        the calculated check points in object space
    %       baseVector:             the base vector in object space
    %       c:                      the focal length in mm
    %       imagePointPrecision:    the image point precision in m
    % Output
    %       differences:    the differences in the real and calculated
    %                       check points
    %       accuracies:     the calculated accuracies of the coordinates

    % Computing the differences
    differences = realCheckPoints - calcCheckPoints;

    H = mean(objectCoordsPC(3, :)); % Computing the flying height from the perspective centers
    h = mean(controlPoints(:, 3));  % Computing the datum height from the control points

    S = (H - h) / (c * 10^-3);  % Computing the scale of the image

    baseMagnitude = sqrt(baseVector(1, 1)^2 + baseVector(2, 1)^2 + baseVector(3, 1)^2); % Computing the magnitude of the base vector
    BH = baseMagnitude / (H - h);                                                       % Computing the base to height ratio

    accuracyX = S * imagePointPrecision;                    % Computing the accuracy in the x-coordinate
    accuracyY = accuracyX;                                  % Computing the accuracy in the y-coordinate
    accuracyZ = ((sqrt(2) * S) / BH) * imagePointPrecision; % Computing the accuracy in the z-coordinate

    accuracies = [accuracyX; accuracyY; accuracyZ]; % Outputting the accuracies into a vector
end
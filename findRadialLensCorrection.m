function [radCorrectionX, radCorrectionY] = findRadialLensCorrection(xBar,yBar)
    % Summary: the findRadialLensCorrection function returns the
    %          corrections to the radial lens distortion in the x and y 
    %          coordinates
    %
    % Input
    %       xBar: reduced x-coordinates to the principle point
    %       yBar: reduced y-coordinates to the principle point
    % Output
    %       radCoorectionX: the corrections to the radial lens distortion
    %                       in the x-coordinate
    %       radCoorectionY: the corrections to the radial lens distortion
    %                       in the y-coordinate

    % Systematic Radial Distortions
    k1 = 0.8878E-4;
    k2 = -0.1528E-7;
    k3 = 0.5256E-12;

    r = sqrt(xBar.^2 + yBar.^2);    % Computing the radial distance

    % Computing the corrections to the radial distortion in the x and y 
    % coordinates 
    radCorrectionX = -xBar .* ((k1 .* r.^2) + (k2 .* r.^4) + (k3 .* r.^6));
    radCorrectionY = -yBar .* ((k1 .* r.^2) + (k2 .* r.^4) + (k3 .* r.^6));
end
function [atmoRefrCorrectionX, atmoRefrCorrectionY] = findAtmosphericRefractionCorrection(xBar, yBar)
    % Summary: the findAtmosphericRefractionCorrection function returns the
    %          corrections to the atmospheric refraction in the x and y
    %          coordinates
    %
    % Input
    %       xBar: reduced x-coordinates to the principle point
    %       yBar: reduced y-coordinates to the principle point
    % Output
    %       atmoRefrCorrectionX: the corrections to the atmospheric
    %                            refraction in the x-coordinate
    %       atmoRefrCorrectionY: the corrections to the atmospheric
    %                            refraction in the y-coordinate
    
    c = 153.358;    % Calibrated Focal Length [mm]

    % Computing the average elevation of the 
    h = mean([1090.96, 1086.43, ...
              1090.50, 1090.65, ...
              1091.55, 1090.82, ...
              1083.49, 1092.00]) / 1000;
    Hh = 751.5 / 1000;  % The flying height of the plane (from Lab 1)
    H = Hh + h;         % The flying height of the plane from the elevation

    r = sqrt(xBar.^2 + yBar.^2);    % Computing the radial distance

    %  Computing the ARDC Model 
    K = ((2410 * H) / (H^2 - 6 * H + 250)) - ((2410 * h) / (h^2 - 6 * h + 250)) * (h / H);

    % Computing the corrections to the atmospheric refraction in the x and
    % y coordinates
    atmoRefrCorrectionX = -xBar .* K .* (1 + (r.^2 ./ c^2)) .* 10E-6;
    atmoRefrCorrectionY = -yBar .* K .* (1 + (r.^2 ./ c^2)) .* 10E-6;
end
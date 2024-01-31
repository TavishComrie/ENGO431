function [xCorr,yCorr] = DecentringLensDistortion(xBar,yBar)
%DecentringLensDistortion to find the decentring lens distortion
% xBar and yBar are provided in fiducial system. Units of mm. May be a vector
% xCorr and yCorr are the corrections in mm for decentring error


    %From the calibration certificate
    p1 = 0.1346e-6;
    p2 = 0.1224e-7;

    %Radial distance
    r = sqrt(xBar.^2+yBar.^2);

    %The correction
    xCorr = -1.*(p1(r.*r+2.*xBar.*xBar)+2.*p2.*xBar.*yBar);
    yCorr = -1.*(p1(r.*r+2.*yBar.*yBar)+2.*p2.*xBar.*yBar);
end
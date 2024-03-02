function [xhat, pY] = performSpaceIntersection(data, dataPrime, x, c)
    bx = 92.000;
    by = x(1, 1);
    bz = x(2, 1);

    for i = 1:size(dataPrime, 1)
        xL = data(i, 1);
        yL = data(i, 2);
    
        xR = dataPrime(i, 1);
        yR = dataPrime(i, 2);
        zR = dataPrime(i, 3);
    
        % Calculating the 
        lambda = (bx * zR - bz * xR) / (xL * zR + c * xR);
        mu = (-bx * c - bz * xL) / (xL * zR + c * xR);

        XL = lambda * xL;
        YL = lambda * yL;
        ZL = -lambda * c;

        XR = mu * xR + bx;
        YR = mu * yR + by;
        ZR = mu * zR + bz;

        X = XL;
        Y = (YL + YR) / 2;
        Z = ZL;

        xhat(i, :) = [X, Y, Z];
        pY = YR - YL;
    end
end
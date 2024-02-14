function [xhat, vhat] = performLeastSquaresAdjustment(data)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    bx = 92.000; 
    xhat = 0;
    vhat = 0;


    function A = findDesignMatrixA(data, dataPrime, xo, C, bx)        
        for i = 1:size(A, 1)
            xL = data(i, 1);
            yL = data(i, 2);

            xR = dataPrime(i, 1);
            yR = dataPrime(i, 2);
            zR = dataPrime(i, 3);

            by = xo(i, 1);
            bz = xo(i, 2);
            omega = xo(i, 3);
            phi = xo(i, 4);

            a = -yR * sin(omega) + zR * cos(omega);
            b = xR * sin(omega);
            c = -xR * cos(omega);
            d = -yR * cos(omega) * cos(phi) - zR * sin(omega) * cos(phi);
            e = xR * cos(omega) * cos(phi) - zR * sin(phi);
            f = xR * sin(omega) * cos(phi) + yR * sin(phi);

            deltaBy = det([0, 1, 0; ...
                           xL, yL, -C; ...
                           xR, yR, zR]);

            deltaBz = det([0, 0, 1; ...
                           xl, yL, -C; ...
                           xR, yR, zR]);

            deltaW = det([bx, by, bz; ...
                          xL, yL, -C; ...
                          0, -zR, yR]);

            deltaPhi = det([bx, by, bz; ...
                             xL, yL, -C; ...
                             a, b, c]);

            deltaKappa = det([bx, by, bz; ...
                              xL, yL, -C; ...
                              d, e, f]);

            A(i, :) = [deltaBy, deltaBz, deltaW, deltaPhi, deltaKappa];
        end         
    end
end
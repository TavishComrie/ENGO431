function [xhat, vhat] = performLeastSquaresAdjustment(data)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    bx = 92.000; 
    xhat = 0;
    vhat = 0;


    function A = findDesignMatrixA(data, dataPrime, xo)
        
        xL = data()
        deltaBy = [0, 1, 0;
                   ]
    end
end
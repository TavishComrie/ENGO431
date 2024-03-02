function xhat = performSpaceIntersection(data, dataPrime, x)
    bx = 92.000;
    by = x(1, 1);
    bz = x(2, 1);
    w = x(3, 1);
    phi = x(4, 1);

    for i = 1:size(dataPrime, 1)
        xL = data(i, 1);
        yL = data(i, 2);
    
        xR = dataPrime(i, 1);
        yR = dataPrime(i, 2);
        zR = dataPrime(i, 3);
    
        lambda = (bx * zR - bz * xR) / (xL * zR + c * xR);

        mu = (-bx * c - bz * xL) / (xL * zR + c * xR);
    end
end
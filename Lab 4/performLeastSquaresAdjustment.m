function [xhat, residuals,Rx, dataprime] = performLeastSquaresAdjustment(data, c)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %In order 7x1, omega,phi,kappa,tx,ty,tz,scale
    xhat(7,1) = zeros;
    %Straigh level so omega,phi stay at 0
    objectbearing = atan2(data(1,5)-data(2,5),data(1,6),data(2,6));
    modelbearing = atan2(data(1,2)-data(2,2),data(1,3),data(2,3));
    objectdistance = sqrt((data(1,5)-data(2,5))^2+(data(1,6)-data(2,6))^2);
    modeldistance = sqrt((data(1,2)-data(2,2))^2+(data(1,3)-data(2,3))^2);



    xhat(3,1) = objectbearing - modelbearing;
    xhat(7,1) = objectdistance / modeldistance;

    robject = [data(1,5);data(1,6);data(1,7)];
    rmodel = [data(1,2);data(1,3);data(1,4)];
    t0 = robject - xhat(7,1) * M_transformation_Matrix(xhat) * rmodel;

    xhat(4,1) = t0(1,1);
    xhat(5,1) = t0(2,1);
    xhat(6,1) = t0(3,1);
    
    threshold = [0.0001;0.0001;0.0001;0.001;0.001;0.001;0.001];
    
    counter = 0;

    notConverged = true;
    while notConverged
        M = M_transformation_Matrix(xhat);

        A = findDesignMatrixA(data, dataprime, xhat, c, bx);
        
        w = createMisclosure(xhat,data,M);

        N = transpose(A) * A;
        u = transpose(A) * w;

        delta = -1 * inv(N) * u;

        xhat = xhat + delta;

        check = delta > threshold;
        notConverged = ismember(1,check);

        counter = counter + 1;

    end

    M = M_transformation_Matrix(xhat);

    w = createMisclosure(xhat,data,M);
    A = findDesignMatrixA(data, dataprime, xhat, c, bx);
    residuals = A * delta + w;

    aPost = transpose(residuals) * residuals;
    Cx = aPost * inv(N);

    Rx = corrcov(Cx);
end

function A = findDesignMatrixA(data, dataPrime, xo, C, bx)        
    for i = 1:6
        xL = data(i, 1);
        yL = data(i, 2);

        xR = dataPrime(i, 1);
        yR = dataPrime(i, 2);
        zR = dataPrime(i, 3);

        by = xo(1, 1);
        bz = xo(2, 1);
        omega = xo(3, 1);
        phi = xo(4, 1);

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
                       xL, yL, -C; ...
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


function w = createMisclosure(x0,data,M)
    points = size(data,1);
    w(points*3,1)=zeros;
    for i = 1:points 
        robject = [data(i,5);data(i,6);data(i,7)];
        rmodel = [data(i,2);data(i,3);data(i,4)];
        t = [x0(4,1);x0(5,1);x0(6,1)];
        scale = x0(7,1);
    
        pointmisclosure = scale * M * rmodel + t - robject;
    
        w(3i-2,1) = pointmisclosure(1,1);
        w(3i-1,1) = pointmisclosure(2,1);
        w(3i,1) = pointmisclosure(3,1);
    end
end


function M = M_transformation_Matrix(Xnot)
    w = Xnot(3,1);
    phi = Xnot(4,1);
    kappa = Xnot(5,1);
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end
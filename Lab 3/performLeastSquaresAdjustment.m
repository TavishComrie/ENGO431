function [xhat, residuals,Rx, dataprime] = performLeastSquaresAdjustment(data)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %In order 5x1: by,bz,w,theta,kappa
    xhat(5,1) = zeros;
    c = 153.358; %mm
    bx = 92.000; 


       
    threshold = [0.001;0.001;0.0001;0.0001;0.0001;];


    xprime = zeros(6,1);
    yprime = zeros(6,1);
    zprime = zeros(6,1);
    
    counter = 0;

    notConverged = true;
    while notConverged
        M = M_transformation_Matrix(xhat);

        %KRISH TODO (set the input vars)   
        for i = 1:6
            mvectorXYZ = [data(i,3), data(i,4), -c];
            [xprime(i),yprime(i),zprime(i)] = transformation(M,mvectorXYZ);
        end

        dataprime = [xprime, yprime, zprime];
        %dataprime(:,1) = dataprime(:,1) + bx;
        %dataprime(:,2) = dataprime(:,2) + xhat(1,1);
        %dataprime(:,3) = dataprime(:,3) + xhat(2,1);


        %RAYMOND TODO (set the input vars)
        A = findDesignMatrixA(data, dataprime, xhat, c, bx)

        
        w = createMisclosure(xhat,data,dataprime,c,bx)

        N = transpose(A) * A;
        u = transpose(A) * w;

        delta = -1 * inv(N) * u;

        xhat = xhat + delta

        check = delta > threshold;
        notConverged = ismember(1,check);

        counter = counter + 1;

    end

    counter


    M = M_transformation_Matrix(xhat);

    %KRISH TODO (set the input vars)   
    for i = 1:6
        mvectorXYZ = [data(i,3), data(i,4), -c];
        [xprime(i),yprime(i),zprime(i)] = transformation(M,mvectorXYZ);
    end

    dataprime = [xprime, yprime, zprime];


    w = createMisclosure(xhat,data,dataprime,c,bx)
    A = findDesignMatrixA(data, dataprime, xhat, c, bx)
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


function w = createMisclosure(x0,data,dataprime,c,bx)
    w(6,1)=zeros;
    for i = 1:6
        misclosure(3,3) = zeros;

        misclosure(1,1) = bx;
        misclosure(1,2:3) = x0(1:2,1);
        misclosure(2,1:2) = data(i,1:2);
        misclosure(2,3) = -1 * c;
        misclosure(3,:) = dataprime(i,:);
        
        w(i,1) = det(misclosure);

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


function [xprime,yprime,zprime] = transformation(M, vector_x_y_z)
    vector_intermediate = transpose(M) * transpose(vector_x_y_z);
    xprime = vector_intermediate(1,1);
    yprime = vector_intermediate(2,1);
    zprime = vector_intermediate(3,1);
end
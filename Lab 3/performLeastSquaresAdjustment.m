function [xhat, residuals,Rx] = performLeastSquaresAdjustment(data)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %In order 5x1: by,bz,w,theta,kappa
    xhat(5,1) = zeros;
    c = 153.358; %mm
    bx = 92.000; 


       
    threshold = [0.001;0.001;0.001;0.001;0.001;];


    xprime = zeros(6,1);
    yprime = zeros(6,1);
    zprime = zeros(6,1);


    notConverged = true;
    while notConverged

        %KRISH TODO (set the input vars)   
        for i = 1:size(data,1)
            mvectorXYZ = [data(i,3), data(i,4), c];
            M = M_transformation_Matrix(xhat);
            [xprime(i),yprime(i),zprime(i)] = transformation(M,mvectorXYZ);
        end

        dataprime = [xprime, yprime, zprime]

        %RAYMOND TODO (set the input vars)
        A = findDesignMatrixA(data, dataprime, xhat, c, bx);

        
        w = createMisclosure(xhat,data,dataprime,c,bx);

        N = transpose(A) * A;
        u = transpose(A) * w;

        delta = -1 * inv(N) * u;

        xhat = xhat + delta;

        check = delta > threshold;
        notConverged = ismember(1,check);
    end

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
    w = Xnot(1,1);
    phi = Xnot(2,1);
    kappa = Xnot(3,1);
    M = [cosd(phi)*cosd(kappa), cosd(w)*sind(kappa)+sind(w)*sind(phi)*cosd(kappa), sind(w)*sind(kappa)-cosd(w)*sind(phi)*cosd(kappa);
        -cosd(phi)*sind(kappa), cosd(w)*cosd(kappa)-sind(w)*sind(phi)*sind(kappa), sind(w)*cosd(kappa)+cosd(w)*sind(phi)*sind(kappa);
        sind(phi), -sind(w)*cosd(phi), cosd(w)*cosd(phi)];
end


function [xprime,yprime,zprime] = transformation(M, vector_x_y_z)
    vector_intermediate = transpose(M)* transpose(vector_x_y_z);
    xprime = vector_intermediate(1,1);
    yprime = vector_intermediate(2,1);
    zprime = vector_intermediate(3,1);
end
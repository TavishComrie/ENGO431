function [xhat, residuals,Rx, dataprime] = performLeastSquaresAdjustment(data, c)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %In order 7x1, omega,phi,kappa,tx,ty,tz
    xhat(7,1) = zeros;
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
        A = findDesignMatrixA(data, dataprime, xhat, c, bx);

        
        w = createMisclosure(xhat,data,dataprime,c,bx);

        N = transpose(A) * A;
        u = transpose(A) * w;

        delta = -1 * inv(N) * u;

        xhat = xhat + delta;

        check = delta > threshold;
        notConverged = ismember(1,check);

        counter = counter + 1;

    end

    M = M_transformation_Matrix(xhat);

    %KRISH TODO (set the input vars)   
    for i = 1:6
        mvectorXYZ = [data(i,3), data(i,4), -c];
        [xprime(i),yprime(i),zprime(i)] = transformation(M,mvectorXYZ);
    end

    dataprime = [xprime, yprime, zprime];


    w = createMisclosure(xhat,data,dataprime,c,bx);
    A = findDesignMatrixA(data, dataprime, xhat, c, bx);
    residuals = A * delta + w;

    aPost = transpose(residuals) * residuals;
    Cx = aPost * inv(N);

    Rx = corrcov(Cx);
end

function A = findDesignMatrixA(data, dataPrime, xo, C, bx)        
    for i = 1:6
        Xm = coords(i,1);
        Ym = coords(i,2);
        Zm = coords(i,3);
 

        omega = xo(1, 1);
        phi = xo(2, 1);
        kappa = xo(3,1);
        tx = xo(4,1);
        ty = xo(5,1);
        yz = xo(6,1);
        scale = xo(7,1);

        %X derivatives
        deltaomegax = ((Ym * scale) * (-sin(omega) * sin(kappa) + cos(omega) * sin(phi) * cos(kappa))) + ((scale*Zm) * (cos(omega) * sin(kappa) + sin(omega) * sin(phi) * cos(kappa)));
        deltaphix = -scale * Xm * sin(phi) * cos(kappa) + scale * Ym * sin(omega) * cos(phi) * cos(kappa)- scale * Zm *cos(omega) * cos(phi) * cos(kappa);
        deltakx = (-scale * Xm * cos(phi) * sin(kappa) ) + ((scale * Ym ) * (cos(omega)*cos(kappa) - sin(omega) * sin(phi) * sin(kappa))) + ((scale * Zm) * (sin(omega)*cos(kappa) + cos(omega)*sin(phi)*sin(kappa)));
        deltatxx = 1;
        deltatyx = 0;
        deltatzx = 0;
        deltascalex = Xm * Mmatrix(1,1) + Ym * Mmatrix(1,2) * Zm * Mmatrix(1,3);

        %derivates Y
        deltaomegay = ((Ym * scale) * (-sin(omega) * cos(kappa) - cos(omega) * sin(phi) * sin(kappa))) + ((scale*Zm) * (cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa)));
        deltaphiy = scale * Xm * sin(phi) * sin(kappa) - scale * Ym * sin(omega) * cos(phi) * sin(kappa) + scale * Zm *cos(omega) * cos(phi) * sin(kappa);
        deltaky = (-scale * Xm * cos(phi) * cos(kappa)) + ((scale * Ym ) * (-cos(omega)*sin(kappa) - sin(omega) * sin(phi) * cos(kappa))) + ((scale * Zm) * (-sin(omega)*sin(kappa) + cos(omega)*sin(phi)*cos(kappa)));
        deltatxy = 0;
        deltatyy = 1;
        deltatzy = 0;
        deltascaley = Xm * Mmatrix(2,1) + Ym * Mmatrix(2,2) * Zm * Mmatrix(2,3);


        %derivates Z
        deltaomegaz = -scale * Ym * cos(omega) * cos(phi) - scale * Zm * sin(omega) * cos(phi);
        deltaphiz = scale * Xm * cos(phi) + scale * Ym * sin(omega) * sin(phi) - scale * Zm *cos(omega) * sin(phi);
        deltakz = 0;
        deltatxz = 0;
        deltatyz = 0;
        deltatzz = 1;
        deltascalez = Xm * Mmatrix(3,1) + Ym * Mmatrix(3,2) * Zm * Mmatrix(3,3);
        
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
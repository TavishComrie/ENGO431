function [xhat, residuals,Rx,M,t,scale,RValues] = performSinglePhotoResection(imageData,objectData,rmse,c,S,rmax)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here    
    %In order 6x1, Xc,Yc,Zc,omega,phi,kappa
    xhat(6,1) = zeros;
    CL = (rmse^2)*eye(size(data,1)*3);
    %initialize the Cl matrix with the precision
    %Make the weight matrix for the adjustment
    P = inv(CL);

    [a,b,dx,dy] = similarityTransform(imageData,objectData);
    Zave = mean(objectData(:,4));


    xhat(1,1) = dx;
    xhat(2,1) = dy;
    xhat(3,1) = c/1000*sqrt(a*a+b*b)+Zave;
    xhat(4,1) = 0;
    xhat(5,1) = 0;
    xhat(6,1) = atan2(b,a);

    %initialize threhold


    xyThreshold = S * rmse / 10;
    omegaphiThreshold = rmse / (10*c);
    kappaThreshold = rmse / (10*rmax);

    angleThreshold = min(omegaphiThreshold,kappaThreshold);

    threshold = [xyThreshold;xyThreshold;0.0001;angleThreshold;angleThreshold;angleThreshold];
  
   
    counter = 0;

    %Run least squares adjustment with defined paramters and functions for
    %all LSA parameters
    notConverged = true;
    while notConverged
        %Find A
        M = M_transformation_Matrix(xhat)
        %Find M

        A = findDesignMatrixA(data, xhat, M)
        %Find A

        w = createMisclosure(xhat,data,M)
        %Find W

        N = transpose(A) * P * A;
        u = transpose(A) * P * w;
        %Find N and U

        cond(N)
        %Find the condition on N

        delta = -1 * (inv(N) * u);
        %Find Delta (corrections for unknowns)

        xhat = xhat + delta
        %Find Xhat (Corrected unknown parameters)

        %Check for convergence against vector of thresholds
        check = delta > threshold;
        notConverged = ismember(1,check);

        counter = counter + 1;
    end
    %post adjustment procedure
    residuals = A * delta + w;
    M = M_transformation_Matrix(xhat);
    t = [xhat(4,1);xhat(5,1);xhat(6,1)];
    scale = xhat(7,1);
    w = createMisclosure(xhat,data,M);
    A = findDesignMatrixA(data,xhat,M);
    aPost = transpose(residuals) *P* residuals / (size(data,1)*3-7);
    %determine correlation and redundancy
    Cx = aPost * inv(N);
    Rx = corrcov(Cx);
    R = eye(size(data,1)*3) - A * inv(A'*P*A) * A' * P;
    RValues = diag(R);
end

function A = findDesignMatrixA(data, xo, Mmatrix) 
%development of the A matrix
    A(size(data,1)*3,7) = zeros;
    for i = 1:size(data,1)
        %initalize model space coords
        Xm = data(i,2);
        Ym = data(i,3);
        Zm = data(i,4);
 

        %extract unknown parameters and store
        omega = xo(1, 1);
        phi = xo(2, 1);
        kappa = xo(3,1);
        tx = xo(4,1);
        ty = xo(5,1);
        yz = xo(6,1);
        scale = xo(7,1);

        %derivative X
        deltaomegax = ((Ym * scale) * ((-sin(omega) * sin(kappa)) + (cos(omega) * sin(phi) * cos(kappa)))) + ((scale*Zm) * ((cos(omega) * sin(kappa)) + sin(omega) * sin(phi) * cos(kappa)));
        deltaphix = (-scale * Xm * sin(phi) * cos(kappa)) + (scale * Ym * sin(omega) * cos(phi) * cos(kappa))- (scale * Zm *cos(omega) * cos(phi) * cos(kappa));
        deltakx = (-scale * Xm * cos(phi) * sin(kappa) ) + ((scale * Ym ) * (cos(omega)*cos(kappa) - sin(omega) * sin(phi) * sin(kappa))) + ((scale * Zm) * (sin(omega)*cos(kappa) + cos(omega)*sin(phi)*sin(kappa)));
        deltatxx = 1;
        deltatyx = 0;
        deltatzx = 0;
        deltascalex = Xm * Mmatrix(1,1) + Ym * Mmatrix(1,2) + Zm * Mmatrix(1,3);

        %derivates Y
        deltaomegay = ((Ym * scale) * (-sin(omega) * cos(kappa) - cos(omega) * sin(phi) * sin(kappa))) + ((scale*Zm) * (cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa)));
        deltaphiy = scale * Xm * sin(phi) * sin(kappa) - scale * Ym * sin(omega) * cos(phi) * sin(kappa) + scale * Zm *cos(omega) * cos(phi) * sin(kappa);
        deltaky = (-scale * Xm * cos(phi) * cos(kappa)) + ((scale * Ym ) * (-cos(omega)*sin(kappa) - sin(omega) * sin(phi) * cos(kappa))) + ((scale * Zm) * (-sin(omega)*sin(kappa) + cos(omega)*sin(phi)*cos(kappa)));
        deltatxy = 0;
        deltatyy = 1;
        deltatzy = 0;
        deltascaley = Xm * Mmatrix(2,1) + Ym * Mmatrix(2,2) + Zm * Mmatrix(2,3);


        %derivates Z
        deltaomegaz = -scale * Ym * cos(omega) * cos(phi) - scale * Zm * sin(omega) * cos(phi);
        deltaphiz = scale * Xm * cos(phi) + scale * Ym * sin(omega) * sin(phi) - scale * Zm *cos(omega) * sin(phi);
        deltakz = 0;
        deltatxz = 0;
        deltatyz = 0;
        deltatzz = 1;
        deltascalez = Xm * Mmatrix(3,1) + Ym * Mmatrix(3,2) + Zm * Mmatrix(3,3);
  
    %Store all deriviates into A matrix to return A
        A(3*i-2, :) = [deltaomegax, deltaphix, deltakx, deltatxx, deltatyx, deltatzx, deltascalex];
        A(3*i-1, :) = [deltaomegay, deltaphiy, deltaky, deltatxy, deltatyy, deltatzy, deltascaley];
        A(3*i, :) = [deltaomegaz, deltaphiz, deltakz, deltatxz, deltatyz, deltatzz, deltascalez];

    end      
end


function w = createMisclosure(x0,data,M)
    %Determine misclosure
    points = size(data,1);
    w(points*3,1)=zeros;
    for i = 1:points 
        %loop to find point misclosure using model, object, translation and
        %scale paramters
        robject = [data(i,5);data(i,6);data(i,7)];
        rmodel = [data(i,2);data(i,3);data(i,4)];
        t = [x0(4,1);x0(5,1);x0(6,1)];
        scale = x0(7,1);
        %point misclosure equation
        pointmisclosure = scale * M * rmodel + t - robject;
        %store and return misclosure values 
        w(3*i-2,1) = pointmisclosure(1,1);
        w(3*i-1,1) = pointmisclosure(2,1);
        w(3*i,1) = pointmisclosure(3,1);
    end
end


function M = M_transformation_Matrix(Xnot)
    %M matrix developed for transforamtion
    %Xhat parameters extracted to be used 
    w = Xnot(1,1);
    phi = Xnot(2,1);
    kappa = Xnot(3,1);
    %initalize M matrix in radians
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end
function [xhat, residuals,Rx,RValues] = performSinglePhotoResection(imageData,objectData,rmse,c,S,rmax)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here    
    %In order 6x1, Xc,Yc,Zc,omega,phi,kappa
    xhat(6,1) = zeros;
    CL = (rmse^2)*eye(size(imageData,1)*2);
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

    xhat

    xhat = [6338.6;3984.6;1453.1;0;0;-18.854*pi/180]

    %initialize threhold


    xyThreshold = S * rmse / 10 / 1000;
    omegaphiThreshold = rmse / (10*c);
    kappaThreshold = rmse / (10*rmax);

    angleThreshold = min(omegaphiThreshold,kappaThreshold);

    threshold = [xyThreshold;xyThreshold;xyThreshold;angleThreshold;angleThreshold;angleThreshold]
  
   
    counter = 0;

    %Run least squares adjustment with defined paramters and functions for
    %all LSA parameters
    notConverged = true;
    while notConverged
        %Find A
        M = M_transformation_Matrix(xhat);
        %Find M

        [A,w] = findDesignMatrixAandW(imageData,objectData,xhat,M,c)

        N = transpose(A) * P * A;
        u = transpose(A) * P * w;
        %Find N and U

        %Find the condition on N

        delta = -1 * (inv(N) * u);
        %Find Delta (corrections for unknowns)

        xhat = xhat + delta;
        %Find Xhat (Corrected unknown parameters)

        %Check for convergence against vector of thresholds
        check = delta > threshold;
        notConverged = ismember(1,check);

        counter = counter + 1;
    end
    %post adjustment procedure
    residuals = A * delta + w;
    aPost = transpose(residuals) *P* residuals / (size(imageData,1)*2-6);
    %determine correlation and redundancy
    Cx = aPost * inv(N);
    Rx = corrcov(Cx);
    R = eye(size(imageData,1)*2) - A * inv(A'*P*A) * A' * P;
    RValues = diag(R);
end

function [A,w] = findDesignMatrixAandW(imageData,objectData,x,M,c) 
%development of the A matrix
    p = size(imageData,1);
    n = 2*p;

    A(n,6) = zeros;
    w(n,1) = zeros;

    Xc = x(1,1);
    Yc = x(2,1);
    Zc = x(3,1);

    phi = asin(M(3,1));
    k = atan2(-M(2,1),M(1,1));
    omega = atan2(-M(3,2),M(3,3));

    for i = 1:p
        Xi = objectData(i,2);
        Yi = objectData(i,3);
        Zi = objectData(i,4);

        xi = imageData(i,2);
        yi = imageData(i,3);


        U = M(1,1)*(Xi-Xc)+M(1,2)*(Yi-Yc)+M(1,3)*(Zi-Zc);
        V = M(2,1)*(Xi-Xc)+M(2,2)*(Yi-Yc)+M(2,3)*(Zi-Zc);
        W = M(3,1)*(Xi-Xc)+M(3,2)*(Yi-Yc)+M(3,3)*(Zi-Zc);

        dxX = -c*(M(3,1)*U-M(1,1)*W)/(W*W);
        dxY = -c*(M(3,2)*U-M(1,2)*W)/(W*W);
        dxZ = -c*(M(3,3)*U-M(1,3)*W)/(W*W);

        dyX = -c*(M(3,1)*V-M(2,1)*W)/(W*W);
        dyY = -c*(M(3,2)*V-M(2,2)*W)/(W*W);
        dyZ = -c*(M(3,3)*V-M(2,3)*W)/(W*W);

        dxw = (-c/(W*W))*((Yi-Yc)*(U*M(3,3)-W*M(1,3))-(Zi-Zc)*(U*M(3,2)-W*M(1,2)));
        dxk = (-c*V)/W;

        dyw = (-c/(W*W))*((Yi-Yc)*(V*M(3,3)-W*M(2,3))-(Zi-Zc)*(V*M(3,2)-W*M(2,2)));
        dyk = (-c*U)/W;

        
        xterm1 = (Xi-Xc)*(-1*W*sin(phi)*cos(k)-U*cos(phi));
        xterm2 = (Yi-Yc)*(W*sin(omega)*cos(phi)*cos(k)-U*sin(omega)*sin(phi));
        xterm3 = (Zi-Zc)*(-1*W*cos(omega)*cos(phi)*cos(k)+U*cos(omega)*sin(phi));

        yterm1 = (Xi-Xc)*(W*sin(phi)*sin(k)-V*cos(phi));
        yterm2 = (Yi-Yc)*(-1*W*sin(omega)*cos(phi)*sin(k)-V*sin(omega)*sin(phi));
        yterm3 = (Zi-Zc)*(W*cos(omega)*cos(phi)*sin(k)+V*cos(omega)*sin(phi));

        dxo = (-c/(W*W)) * (xterm1+xterm2+xterm3);
        dyo = (-c/(W*W)) * (yterm1+yterm2+yterm3);


        A(2*i-1,1)=dxX;
        A(2*i-1,2)=dxY;
        A(2*i-1,3)=dxZ;
        A(2*i-1,4)=dxw;
        A(2*i-1,5)=dxo;
        A(2*i-1,6)=dxk;

        A(2*i,1)=dyX;
        A(2*i,2)=dyY;
        A(2*i,3)=dyZ;
        A(2*i,4)=dyw;
        A(2*i,5)=dyo;
        A(2*i,6)=dyk;

        fx = -c*U/W;
        fy = -c*V/W;   

        wx = fx - xi;
        wy = fy - yi;

        w(2*i-1,1)=wx;
        w(2*i,1)=wy;

    end
    w
end



function M = M_transformation_Matrix(Xnot)
    %M matrix developed for transforamtion
    %Xhat parameters extracted to be used 
    w = Xnot(4,1);
    phi = Xnot(5,1);
    kappa = Xnot(6,1);
    %initalize M matrix in radians
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end
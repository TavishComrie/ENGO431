function [xhat, residuals,Rx,M,t,scale,RValues] = performCollinearityLSA(data28, EOP, IOP, data27 )
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here    
    %In order 7x1, omega,phi,kappa,tx,ty,tz,scale
    xhat(3,1) = zeros;
    CL = eye(4,4);
    %initialize the Cl matrix with the precision

    %Precision of the x,y coords in image space
    %is taken from Lab 2 report as RMS the values
    precImg = 0.004;

    CL = (1/precImg) * CL;
    %Make the weight matrix for the adjustment
    P = inv(CL);
    
    Xc = EOP(1,:);
    Yc = EOP(2,:);
    Zc = EOP(3,:);

    EOPAngles = EOP(4:6,:);

    
    %initialize threhold
    threshold = [0.0001;0.0001;0.0001;0.001;0.001;0.001;0.001];
    
    xo(size(data28,1),1) = zeros;
    yo(size(data28,1),1) = zeros;
    
    %Find approximate Values

    xVals28 = data28(:,1);
    yVals28 = data28(:,2);
    
    xVals27 = data27(:,1);
    yVals28 = data28(:,1);

    for i = 1:size(data28,1)
        for j = 1:2
            x(

        end
    end
    

   
    counter = 0;

    %Run least squares adjustment with defined paramters and functions for
    %all LSA parameters
    notConverged = true;
    while notConverged
        
        M27 = M_transformation_Matrix(EOPAngles(:,1))
        M28 = M_transformation_Matrix(EOPAngles(:,2))
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


function [A,w] = findDesignMatrixAandW(xhat,imagePoints,ML,MR,c,EOP) 
%development of the A matrix
    n = 4;

    A(4,3) = zeros;
    w(4,1) = zeros;

    Xcl = EOP(1,1);
    Ycl = EOP(2,1);
    Zcl = EOP(3,1);

    Xcr = EOP(1,2); 
    Ycr = EOP(2,2);
    Zcr = EOP(3,2);


    Xi = xhat(1,1); 
    Yi = xhat(2,1);
    Zi = xhat(3,1);

    xiL = imageData(1,2);
    yiL = imageData(1,3);

    xiR = imageData(1,4);
    yiR = imageData(1,5);


    Ul = ML(1,1)*(Xi-Xcl)+ML(1,2)*(Yi-Ycl)+ML(1,3)*(Zi-Zcl);
    Vl = ML(2,1)*(Xi-Xcl)+ML(2,2)*(Yi-Ycl)+ML(2,3)*(Zi-Zcl);
    Wl = ML(3,1)*(Xi-Xcl)+ML(3,2)*(Yi-Ycl)+ML(3,3)*(Zi-Zcl);

    Ur = MR(1,1)*(Xi-Xcr)+MR(1,2)*(Yi-Ycr)+MR(1,3)*(Zi-Zcr);
    Vr = MR(2,1)*(Xi-Xcr)+MR(2,2)*(Yi-Ycr)+MR(2,3)*(Zi-Zcr);
    Wr = MR(3,1)*(Xi-Xcr)+MR(3,2)*(Yi-Ycr)+MR(3,3)*(Zi-Zcr);

        A(1,1) = c*(ML(3,1)*U-M(1,1)*W)/(W*W);
        A(1,2) = c*(M1(3,2)*U-M(1,2)*W)/(W*W);
        A(1,3) = c*(M1(3,3)*U-M(1,3)*W)/(W*W);

        dyX = c*(M(3,1)*V-M(2,1)*W)/(W*W);
        dyY = c*(M(3,2)*V-M(2,2)*W)/(W*W);
        dyZ = c*(M(3,3)*V-M(2,3)*W)/(W*W);


        A(2*i-1,1)=dxX;
        A(2*i-1,2)=dxY;
        A(2*i-1,3)=dxZ;

        A(2*i,1)=dyX;
        A(2*i,2)=dyY;
        A(2*i,3)=dyZ;

        fx = -c*U/W;
        fy = -c*V/W;   

        wx = fx - xi;
        wy = fy - yi;

        w(2*i-1,1)=wx;
        w(2*i,1)=wy;

    end
    w
end


function M = M_transformation_Matrix(EOP)
    %M matrix developed for transforamtion
    %Xhat parameters extracted to be used 
    w = EOP(4,1);
    phi = EOP(5,1);
    kappa = EOP(6,1);
    %initalize M matrix in radians
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end
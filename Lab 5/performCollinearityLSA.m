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
        dyk = (c*U)/W;

        
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
    w = Xnot(1,1);
    phi = Xnot(2,1);
    kappa = Xnot(3,1);
    %initalize M matrix in radians
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end
function [xhat, residuals,Rx] = performLeastSquaresAdjustment(data)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %In order 5x1: by,bz,w,theta,kappa
    xhat(5,1) = zeros;
    c = 153.358; %mm

    threshold = [0.001;0.001;0.001;0.001;0.001;];

    notConverged = false;
    while notConverged

        %KRISH TODO (set the input vars)
        M = M_transformation_matrix();
        dataprime = transformation();

        %RAYMOND TODO (set the input vars)
        A = findDesignMatrix();

        
        w = createMisclosure(x0,data,dataprime,c,bx);

        N = transpose(A) * A;
        u = transpose(A) * w;

        delta = inv(N) * u;

        xhat = xhat + delta;

        check = delta > threshold;
        notConverged = ismember(1,check);
    end

    residuals = A * delta + w;

    aPost = transpose(residuals) * residuals;
    Cx = aPost * inv(N);

    Rx = corrcov(Cx);
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
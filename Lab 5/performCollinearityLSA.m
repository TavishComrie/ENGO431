function [xhat, redundancyNumbers] = performCollinearityLSA(c,data28, EOP, data27 )
    % Summary:  collinearityLSA transforms
    %
    % Input
    %       c:          the focal length 
    %       data28:     image point coords of image 28
    %       data27:     image point coords of image 27
    %       EOP:        EOPs of each image
    % Output
    %       xhat:               the estimated EOPs
    %       redundancyNumbers:  the vector of redundancy values
    CL = eye(4,4);
    %initialize the Cl matrix with the precision

    %Precision of the x,y coords in image space
    %is taken from Lab 2 report as RMS the values
    precImg = 0.004;
    
    CL = (precImg) * CL;
    %Make the weight matrix for the adjustment
    P = inv(CL);
    
    %Reads the EOP values
    Xc = EOP(1,:);
    
    Yc = EOP(2,:);
    
    Zc = EOP(3,:);

    EOPAngles = EOP(4:6,:);
  
    %initialize threhold
    threshold = [0.0001;0.0001;0.0001];
    

    %Find approximate Values

    xVal28 = data28(:,1);
    yVal28 = data28(:,2);
    
    xVal27 = data27(:,1);
    yVal27 = data27(:,2);


    M27 = M_transformation_Matrix(EOPAngles(:,1));
    M28 = M_transformation_Matrix(EOPAngles(:,2));

    weirdA = make_weirdA(c,M27,M28,xVal27,yVal27,xVal28,yVal28);
    b = make_b(c,M27,M28,xVal27,yVal27,xVal28,yVal28,Xc,Yc,Zc);

    weirdA
    b
   
   
    %Approximate LSA runs
    xhat = inv((transpose(weirdA) * weirdA))*transpose(weirdA)*b;
    xhat
    counter = 0;

    %Approximate values solved for, starting adjustment

    %Run least squares adjustment with defined paramters and functions for
    %all LSA parameters
    notConverged = true;
    while notConverged
        
       
        %Find M

        [A,w] = findDesignMatrixAandW(xhat,data27,data28,M27,M28,c,EOP); 

        N = transpose(A) * P * A;
        u = transpose(A) * P * w;
        %Find N and U

        cond(N);
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
    
    I = eye(size(data27, 2) * 2);
    R = I - (A * inv(A'*P*A) * A' * P);
    redundancyNumbers = diag(R);
end


function [A,w] = findDesignMatrixAandW(xhat,dataL,dataR,ML,MR,c,EOP) 
%development of the A matrix
    n = 4;

    %Sizes A and w
    A(4,3) = zeros;
    w(4,1) = zeros;

    %Reads EOP's
    Xcl = EOP(1,1);
    Ycl = EOP(2,1);
    Zcl = EOP(3,1);

    Xcr = EOP(1,2); 
    Ycr = EOP(2,2);
    Zcr = EOP(3,2);

    %Reads unknowns
    Xi = xhat(1,1); 
    Yi = xhat(2,1);
    Zi = xhat(3,1);

    xiL = dataL(1,1);
    yiL = dataL(1,2);

    xiR = dataR(1,1);
    yiR = dataR(1,2);

    %Applies collinearity equations
    Ul = ML(1,1)*(Xi-Xcl)+ML(1,2)*(Yi-Ycl)+ML(1,3)*(Zi-Zcl);
    Vl = ML(2,1)*(Xi-Xcl)+ML(2,2)*(Yi-Ycl)+ML(2,3)*(Zi-Zcl);
    Wl = ML(3,1)*(Xi-Xcl)+ML(3,2)*(Yi-Ycl)+ML(3,3)*(Zi-Zcl);

    Ur = MR(1,1)*(Xi-Xcr)+MR(1,2)*(Yi-Ycr)+MR(1,3)*(Zi-Zcr);
    Vr = MR(2,1)*(Xi-Xcr)+MR(2,2)*(Yi-Ycr)+MR(2,3)*(Zi-Zcr);
    Wr = MR(3,1)*(Xi-Xcr)+MR(3,2)*(Yi-Ycr)+MR(3,3)*(Zi-Zcr);

    %Applies derivatives
    A(1,1) = c*(ML(3,1)*Ul-ML(1,1)*Wl)/(Wl*Wl);
    A(1,2) = c*(ML(3,2)*Ul-ML(1,2)*Wl)/(Wl*Wl);
    A(1,3) = c*(ML(3,3)*Ul-ML(1,3)*Wl)/(Wl*Wl);

    A(2,1) = c*(ML(3,1)*Vl-ML(2,1)*Wl)/(Wl*Wl);
    A(2,2) = c*(ML(3,2)*Vl-ML(2,2)*Wl)/(Wl*Wl);
    A(2,3) = c*(ML(3,3)*Vl-ML(2,3)*Wl)/(Wl*Wl);

    A(3,1) = c*(MR(3,1)*Ur-MR(1,1)*Wr)/(Wr*Wr);
    A(3,2) = c*(MR(3,2)*Ur-MR(1,2)*Wr)/(Wr*Wr);
    A(3,3) = c*(MR(3,3)*Ur-MR(1,3)*Wr)/(Wr*Wr);

    A(4,1) = c*(MR(3,1)*Vr-MR(2,1)*Wr)/(Wr*Wr);
    A(4,2) = c*(MR(3,2)*Vr-MR(2,2)*Wr)/(Wr*Wr);
    A(4,3) = c*(MR(3,3)*Vr-MR(2,3)*Wr)/(Wr*Wr);


    %Solves for misclosure
    fxL = -c*Ul/Wl;
    fyL = -c*Vl/Wl;   

    fxR = -c*Ur/Wr;
    fyR = -c*Vr/Wr;

    wxL = fxL - xiL;
    wyL = fyL - yiL;

    wxR = fxR - xiR;
    wyR = fyR - yiR;

    w(1, 1)=wxL;
    w(2, 1)=wyL;
    w(3, 1)=wxR;
    w(4, 1)=wyR;

    
end


function M = M_transformation_Matrix(EOP)
    %M matrix developed for transforamtion
    %Xhat parameters extracted to be used 
    w = EOP(1,1);
    phi = EOP(2,1);
    kappa = EOP(3,1);
    %initalize M matrix in radians
    M = [cos(phi)*cos(kappa), cos(w)*sin(kappa)+sin(w)*sin(phi)*cos(kappa), sin(w)*sin(kappa)-cos(w)*sin(phi)*cos(kappa);
        -cos(phi)*sin(kappa), cos(w)*cos(kappa)-sin(w)*sin(phi)*sin(kappa), sin(w)*cos(kappa)+cos(w)*sin(phi)*sin(kappa);
        sin(phi), -sin(w)*cos(phi), cos(w)*cos(phi)];
end


%Functions to find approximate values
function weirdA = make_weirdA(c,mL,mR,xL,yL,xR,yR)


        weirdA(1,1) = (xL(1)*mL(3,1))+(c*mL(1,1));
        weirdA(1,2) = (xL(1)*mL(3,2))+(c*mL(1,2));
        weirdA(1,3) = (xL(1)*mL(3,3))+(c*mL(1,3));

        weirdA(2,1) = (yL(1)*mL(3,1))+(c*mL(2,1));
        weirdA(2,2) = (yL(1)*mL(3,2))+(c*mL(2,2));
        weirdA(2,3) = (yL(1)*mL(3,3))+(c*mL(2,3));

        weirdA(3,1) = (xR(1)*mR(3,1))+(c*mR(1,1));
        weirdA(3,2) = (xR(1)*mR(3,2))+(c*mR(1,2));
        weirdA(3,3) = (xR(1)*mR(3,3))+(c*mR(1,3));

        weirdA(4,1) = (yR(1)*mR(3,1))+(c*mR(2,1));
        weirdA(4,2) = (yR(1)*mR(3,2))+(c*mR(2,2));
        weirdA(4,3) = (yR(1)*mR(3,3))+(c*mR(2,3));
   
end

function b = make_b(c,mL,mR,xL,yL,xR,yR,Xc,Yc,Zc)
        b(4,1) = zeros;

         b(1,1) = (xL(1)*mL(3,1)+(c*mL(1,1)))*Xc(1,1)...
         + ((xL(1)*mL(3,2)+(c*mL(1,2)))*Yc(1,1))...
         + ((xL(1)*mL(3,3))+(c*mL(1,3)))*Zc(1,1);

         b(2,1) = (yL(1)*mL(3,1)+(c*mL(2,1))) *Xc(1,1)...
         + ((yL(1)*mL(3,2)+(c*mL(2,2)))*Yc(1,1))...
         + ((yL(1)*mL(3,3))+(c*mL(2,3)))*Zc(1,1);

         b(3,1) = (xR(1)*mR(3,1)+(c*mR(1,1))) *Xc(1,2)...
         + ((xR(1)*mR(3,2)+(c*mR(1,2))))*Yc(1,2)...
         + ((xR(1)*mR(3,3))+(c*mR(1,3)))*Zc(1,2);

         b(4,1) = (yR(1)*mR(3,1)+(c*mR(2,1))) *Xc(1,2)...
         + ((yR(1)*mR(3,2)+(c*mR(2,2))))*Yc(1,2)...
         + ((yR(1)*mR(3,3))+(c*mR(2,3)))*Zc(1,2);
end




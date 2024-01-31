fids = [1346.25	-19284.75	1345.75	-19286
19161	-1470.25	19174.75	-1483.25
1345	-1471.75	1358.25	-1471.75
19162.5	-19283.75	19162.25	-19297.75
842	-10378	848.75	-10378.25
19665.25	-10376.5	19672.25	-10391.25
10254	-966.25	10267	-973
10254.25	-19788.25	10253.25	-19795.5
];
obsin = [-105.997	-105.995
106.004	106.008
-106	106.009
106.012	-105.995
-112	0.007
112.006	0.007
0.005	112.007
0.002	-111.998
];

obs = normalize(obsin(:,1),obsin(:,2));

fids27 = normalize(fids(:,1),fids(:,2));
fids28 = normalize(fids(:,3),fids(:,4));


A27 = makeA(fids27);
A28 = makeA(fids28);

xhat27 = makeXhat(A27,obs);
xhat28 = makeXhat(A28,obs);

v27 = residuals(A27,xhat27,obs);
v28 = residuals(A28,xhat28,obs);


[theta,Sx,Sy,skew] = calcParams(xhat27(1),xhat27(2),xhat27(4),xhat27(5))











function A = makeA(fids)
    A = zeros(size(fids,1),6);
    for i = 1:2:size(fids,1)
            A(i,1) = fids(i);
            A(i,2) = fids(i+1);
            A(i,3) = 1;
            A(i+1,4) = fids(i);
            A(i+1,5) = fids(i+1);
            A(i+1,6) = 1;
            
    end
end


function xhat = makeXhat(A,obs)
   xhat = ((transpose(A)*A)^-1)*transpose(A)*obs;
end

function v = residuals(A,xhat,obs)
    v = (A*xhat)-obs;
end

function [theta,Sx,Sy,skew] = calcParams(a,b,c,d)
    theta = atan(c/a);
    Sx = sqrt(a^2+c^2);
    Sy = sqrt(b^2+d^2);
    skew = atan((a*b+c*d)/(a*d-b*d));
end

function normal = normalize(col1,col2)
    normal = zeros(size(col1,1)*2,1);
    j=1;
    
    for i = 1:2:size(col1,1)*2
        normal(i) = col1(j);
        normal(i+1) = col2(j);
    
        j= j+1;
    end
end
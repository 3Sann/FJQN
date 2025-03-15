lamda = 10/3600; 
N = 1;
K = 6;
P = 6;

temp = 1;
for k = 1:10
    inTemp =1;
    for i =1:k
        newOutTask = k_equal(i,N);
        TH = OutMVA(K,P,newOutTask);
        inTemp = inTemp * (lamda/TH);
    end
    temp = temp + inTemp;
end
pai(1) = 1/temp;
for k = 1:10
    inTemp=1;
    for i = 1:k
        newOutTask = k_equal(i,N);
        TH = OutMVA(K,P,newOutTask);
        inTemp = inTemp * (lamda/TH);
    end
    pai(k+1) = pai(1)*inTemp;
end
K1 = 0;
u1=0;
for k = 1:10
    K1 = K1 + k * pai(k+1); 
    newOutTask = k_equal(k,N);
    TH = OutMVA(K,P,newOutTask);
    u1 = u1 +  TH * pai(k+1);
end

T = K1/lamda
W = T - 1/u1
Leq = lamda*W

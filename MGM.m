% The set of parameters
% The number of robots
K=6;
% The number of pickers
P=6;
% The number of zones
N=1;
% The arrival rates of order batches
lamda = 60/3600; 
% the service time of the second station
mu2 = 10;

% the computation of B00
B00 = zeros((K+1)*(K+2)/2,(K+1)*(K+2)/2);
for i = 0:K
    if i == 0
        B00(1,1) = -namuta;
        B00(1,2) = namuta;
    else
        C0 = zeros(i+1,i);
        for j = 1:i
            C0(j+1,j) = j/mu2;
        end
        B0 = zeros(i+1,i+1);
        temp = size(B0);
        for j = 1:temp(1)
            if j == 1
                kequal = k_equal(i-j+1,N);
                u = AVMA(K,P,kequal);
                B0(1,1) = - namuta - u;
                B0(1,2) = u;
            elseif j == temp(1)
                B0(j,j) = -(j-1)/mu2 - namuta;
            else
                kequal = k_equal(i-j+1,N);
                 u = AVMA(K,P,kequal);
                B0(j,j) = - namuta - u - (j-1)/mu2;
                B0(j,j+1) = u;
            end
        end
        if i ~= K
            A0 = namuta*eye(i+1);
            Z = [C0 B0 A0];
        else
            Z = [C0 B0];
        end
        temp = size(Z);
        B00((i+1)*i/2+1:(i+1)*i/2+temp(1),(i-1)*i/2+1:(i-1)*i/2+temp(2)) = Z;
    end
end


% The computation of B01
B01 = zeros((K+1)*(K+2)/2,K+1);
temp = size(B01);
B01(temp(1)-K:temp(1),:)= namuta*eye(K+1);

% The computation of B10
B10 = zeros(K+1,(K+1)*(K+2)/2);
C = zeros(K+1,K+1);
for i = 2:K+1
    C(i,i-1) = (i-1)/mu2;
end
temp=size(B10);
B10(:,temp(2)-K:temp(2)) = C;

A = namuta*eye(K+1);

% The computation of B
B = zeros(K+1,K+1);
for i = 0:K
    if i <K
        kequal = k_equal(K-i,N);
        u = AVMA(K,P,kequal);
    end
    if i == 0
        B(i+1,i+1) = - namuta - u;
        B(i+1,i+2) = u;
    elseif i<K
        B(i+1,i+1) = - namuta - u - i/mu2;
        B(i+1,i+2) = u;
    else
        B(i+1,i+1) = - namuta - i/mu2;
    end
end

R=zeros(K+1,K+1);
R2=-(A+R^2*C)/B;
while norm(R2-R,'inf') > 0.01
    R=R2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    R2=-(A+R^2*C)/B;
end
R=R2;

I = eye(K+1);
F = inv(I - R);

e_piao = ones((K+1)*(K+2)/2, 1);
e = ones(K+1,1);
firstline = [e_piao B00 B01];
secondline = [(I - R)\e B10 B+R*C]; 
prodMatrix = [firstline;secondline];
rightMatrix = zeros(1, 1 + (K+1)*(K+2)/2 + (K+1));
rightMatrix(1,1) = 1;
pai = rightMatrix/prodMatrix;

pi0 = pai(1:(K+1)*(K+2)/2);
pi1 = pai((K+1)*(K+2)/2+1:(K+1)*(K+2)/2 + (K+1));

alpha = sum(pi0) + sum(pi1/(I - R));

pi0 = pi0/alpha;
pi1 = pi1/alpha;

Leq = pi1*F^2*ones(K+1,1);

LDC=0;
for k=1:K
    for i=(k+1)*k/2+1: (k+1)*k/2+k+1
        LDC= LDC + k*pi0(1,i);
    end
end

LDC=LDC+N*pi1*F*ones(K+1,1);
W1 = Leq/namuta;
W2 = LDC/namuta;
T=(Leq+LDC)/namuta;
pRo = LDC/K;

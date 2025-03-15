% The set of parameters
% The number of resources
N = 6;
% The arrival rates of order batches
lamda = 60/3600; 
mu2 = 10;

% the computation of B00
B00 = zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);
for i = 0:N
    if i == 0
        B00(1,1) = -lamda;
        B00(1,2) = lamda;
    else
        C0 = zeros(i+1,i);
        for j = 1:i
            C0(j+1,j) = j/mu2;
        end
        B0 = zeros(i+1,i+1);
        temp = size(B0);
        for j = 1:temp(1)
            if j == 1
                kequal = k_equal(i-j+1);
                u = AVMA(kequal);
                B0(1,1) = -lamda - u;
                B0(1,2) = u;
            elseif j == temp(1)
                B0(j,j) = -(j-1)/mu2 - lamda;
            else
                kequal = k_equal(i-j+1);
                u = AVMA(kequal);
                B0(j,j) = -lamda - u - (j-1)/mu2;
                B0(j,j+1) = u;
            end
        end
        if i ~= N
            A0 = lamda*eye(i+1);
            Z = [C0 B0 A0];
        else
            Z = [C0 B0];
        end
        temp = size(Z);
        B00((i+1)*i/2+1:(i+1)*i/2+temp(1),(i-1)*i/2+1:(i-1)*i/2+temp(2)) = Z;
    end
end

% The computation of B01
B01 = zeros((N+1)*(N+2)/2,N+1);
temp = size(B01);
B01(temp(1)-N:temp(1),:)= lamda*eye(N+1);

% The computation of B10
B10 = zeros(N+1,(N+1)*(N+2)/2);
C = zeros(N+1,N+1);
for i = 2:N+1
    C(i,i-1) = (i-1)/mu2;
end
temp=size(B10);
B10(:,temp(2)-N:temp(2)) = C;

A = namuta*eye(N+1);

% The computation of B
B = zeros(N+1,N+1);
for i = 0:N
    if i <N
        kequal = k_equal(N-i);
        u = AVMA(kequal);
    end
    if i == 0
        B(i+1,i+1) = - lamda - u;
        B(i+1,i+2) = u;
    elseif i<N
        B(i+1,i+1) = - lamda - u - i/mu2;
        B(i+1,i+2) = u;
    else
        B(i+1,i+1) = - lamda - i/mu2;
    end
end


R=zeros(N+1,N+1);
R2=-(A+R^2*C)/B;
while norm(R2-R,'inf') > 0.01
    R=R2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    R2=-(A+R^2*C)/B;
end
R=R2;

I = eye(N+1);
F = inv(I - R);

e_piao = ones((N+1)*(N+2)/2, 1);
e = ones(N+1,1);
firstline = [e_piao B00 B01];
secondline = [(I - R)\e B10 B+R*C]; 
prodMatrix = [firstline;secondline];
rightMatrix = zeros(1, 1 + (N+1)*(N+2)/2 + (N+1));
rightMatrix(1,1) = 1;
pai = rightMatrix/prodMatrix;


pi0 = pai(1:(N+1)*(N+2)/2);
pi1 = pai((N+1)*(N+2)/2+1:(N+1)*(N+2)/2 + (N+1));

alpha = sum(pi0) + sum(pi1/(I - R));

pi0 = pi0/alpha;
pi1 = pi1/alpha;

Leq = pi1*F^2*ones(N+1,1);

LDC=0;
for k=1:N
    for i=(k+1)*k/2+1: (k+1)*k/2+k+1
        LDC= LDC + k*pi0(1,i);
    end
end

LDC=LDC+N*pi1*F*ones(N+1,1);
W1 = Leq/lamda;
W2 = LDC/lamda;
T=(Leq+LDC)/lamda;
pRo = LDC/N;

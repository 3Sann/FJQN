% Parameters
vp=0.75;vr=3;
n=64;m=30;na=16;
r=n/na;
f=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
wl=1;
wc=3;
C1=(2*wl+wc)/f; 
Cv=vr/vp;
pTperItem = 3;
pickTime = n*pTperItem;

% The zeros array 
prob = zeros(1,m);
probr = zeros(1,m);

for i=1:m
    prob(i)=((m+1-i)^r-(m-i)^r)/m^r;
end
for i=1:m
    probr(i) = (i^r -(i-1)^r)/m^r;
end

% Go back to the buffer
Eb = 0;
for i=0:m-1
    Eb = Eb + probr(i+1)*(i+wl);
end
Eb = Eb/vr;


% Go to the first picking position (robots)
Efr=0;
for i=0:m-1
    Efr = Efr + prob(i+1)*k;
end
Efr = (Efr + 1/2*(na-1)*wc)/vr;

% Go to the first picking postion (pickers)
Efp=0;
Efp2 = 0;
for i=1:m
   Efp = Efp + prob(i) * abs(m-2*i)/vp;
end


% Go back to the depot
Ed=0;
for i =0:m-1
    Ed = Ed + probr(i+1)*i;
end
Ed = (Ed + wl + 1/2*(na-1)*wc)/vr;


% From the first picking position to the last picking positions
temp = 0;
for i=0:m-1
    temp = temp + probr(i+1)*(i+wl);
end

Elp = ((na-1)*wc + (2*na-1)*temp)/vp;
Elr = ((na-1)*wc + (2*na-1)*temp)/vr;


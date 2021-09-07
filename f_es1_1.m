function dy = f_es1_1(t,y,k)

dy=zeros(2,1);
k1 = k(1);
kmeno1 = k(2);
k2 = k(3);
e0 = k(4);
s = y(1);
c = y(2);
dy(1) = kmeno1*c-k1*s*(e0-c);  %s
dy(2) = k1*s*(e0-c)-(kmeno1+k2)*c; %c
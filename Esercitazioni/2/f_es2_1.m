function dy = f_es2_1(t,y,k)

k1 = k(1);
k_meno1 = k(2);
k2 = k(3);
k3 = k(4);
k_meno3 = k(5);
k4 = k(6);
e0 = k(7);


dy = zeros(3,1);

s = y(1);
c1 = y(2);
c2 = y(3);

%%%%%%%%%%%%%%%%%%%%%%%
e = e0 - c1 - c2;
%%%%%%%%%%%%%%%%%%%%
dy(1) = -k1*e*s + k_meno1*c1 - k3*c1*s + k_meno3*c2;
dy(2) = k1*e*s - (k_meno1 + k2)*c1 - k3*c1*s + (k_meno3+k4)*c2;
dy(3) = k3*c1*s - (k_meno3+k4)*c2;

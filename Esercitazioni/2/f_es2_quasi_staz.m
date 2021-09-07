function ds = f_es2_quasi_staz(t,s,k)

k1 = k(1);
k_meno1 = k(2);
k2 = k(3);
k3 = k(4);
k_meno3 = k(5);
k4 = k(6);
e0 = k(7);
K1 = k(8);
K2 = k(9);


ds = 0;
c1 = K2*e0*s./(K1*K2 + K2*s + s.^2);
c2 = e0*s.^2./(K1*K2 + K2*s + s.^2);
e = e0 - c1 - c2;

ds = -k1*e.*s + k_meno1*c1 - k3*c1.*s + k_meno3*c2;

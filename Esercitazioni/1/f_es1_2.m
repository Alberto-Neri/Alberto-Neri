function dyu = f_es1_2(t,yu,k)

dyu = 0;

k1 = k(1);
kmeno1 = k(2);
k2 = k(3);
e0 = k(4);
Km = (kmeno1+k2)/k1;

dyu = - k2*e0*yu/(Km+yu);
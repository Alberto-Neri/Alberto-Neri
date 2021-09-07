function dy = f_es3_1_quasi_staz(t,y,k)

dy = zeros(2,1);

km = (k(2)+k(3))/k(1);
ki = k(5)/k(4);
c1 = k(6)*y(1)/(y(1)+km*(1+y(2)/ki));
c2 = k(6)*y(2)/(y(2)+ki*(1+y(1)/km));
e = k(6) - c1 -c2;
dy(1) = - k(1)*y(1)*e + k(2)*c1; %s
dy(2) = - k(4)*y(2)*e + k(5)*c2; %i
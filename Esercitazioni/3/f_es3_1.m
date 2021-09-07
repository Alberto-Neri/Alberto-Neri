function dy = f_es3_1(t,y,k)

dy = zeros(4,1);
e = k(6) - y(3) - y(4);

dy(1) = - k(1)*y(1)*e + k(2)*y(3); %s
dy(2) = - k(4)*y(2)*e + k(5)*y(4); %i
dy(3) = k(1)*y(1)*e - (k(2)+k(3))*y(3); %c1
dy(4) = k(4)*y(2)*e - k(5)*y(4);  %c2
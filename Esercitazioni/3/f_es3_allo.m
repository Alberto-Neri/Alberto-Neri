function dy = f_es3_allo(t,y,k)

dy = zeros(5,1);
e = k(6) - y(2) - y(3) - y(4);

dy(1) = -e*y(1)*k(1) + k(2)*y(2); % ds
dy(2) = k(1)*y(1)*e - (k(2)+k(3)+k(4)*y(5))*y(2) + k(5)*y(4); %dx 
dy(3) = k(4)*y(5)*e - (k(5)+k(1)*y(1))*y(3) + k(2)*y(4); %dy
dy(4) = k(1)*y(1)*y(3) - (k(2)+k(5))*y(4) + k(4)*y(5)*y(2); % dz
dy(5) = -e*y(5)*k(4) + k(5)*y(3); %di
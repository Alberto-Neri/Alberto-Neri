function dy = f_es3_allo_qs(t,y,k)

dy = zeros(2,1);

k1 = k(1);
km1 = k(2);
k2 = k(3);
k3 = k(4);
km3 = k(5);
e0 = k(6);

s = y(1);
i = y(2);
x = (e0*k1*km3*s*(km1 + km3 + i*k3 + k1*s))/((km3 + i*k3)*(k2*km1 + k2*km3 + km1*km3 + km1^2 + k1^2*s^2 + i*k3*km1 + k1*k2*s + 2*k1*km1*s + k1*km3*s + i*k1*k3*s));
y = (e0*i*k3*(k2*km1 + k2*km3 + km1*km3 + km1^2 + i*k3*km1 + k1*km1*s))/((km3 + i*k3)*(k2*km1 + k2*km3 + km1*km3 + km1^2 + k1^2*s^2 + i*k3*km1 + k1*k2*s + 2*k1*km1*s + k1*km3*s + i*k1*k3*s));
z = (e0*i*k1*k3*s*(k2 + km1 + km3 + i*k3 + k1*s))/((km3 + i*k3)*(k2*km1 + k2*km3 + km1*km3 + km1^2 + k1^2*s^2 + i*k3*km1 + k1*k2*s + 2*k1*km1*s + k1*km3*s + i*k1*k3*s));

e = e0 - x - y - z;

dy(1) = -e*s*k1 + km1*x;
dy(2) = -e*i*k3 + km3*y;
function dy = fhn_es6(t,y,k)

dy = zeros(2,1); % v w
b = k(2);
beta = k(4);
delta = k(5);
c = k(3);
Iapp = k(1);
gamma = k(6);
e = k(7);

if t>0.5  %do impulso solo fino al tempo che voglio 
    Iapp =0; %Iapp
end
% k = Iapp b c beta delta gamma e

dy(1) = b*y(1)*(y(1)-beta)*(delta-y(1))-c*y(2)+Iapp;
dy(2) = e*(y(1)-gamma*y(2));

end


clear all
close all
clc

%% inibizione competitiva

t = [0,100];
y0 = [5,2,0,0];
k = [102,50,1,26,50,1];

[tout_ex yout] = ode15s(@f_es3_1,t,y0,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,4),'InitialStep',1.e-5),k);

V = k(3)*yout(:,3);

figure(1)
plot(yout(:,1),V)
axis([0, 2, 0, inf])
title('V(s) inib competitiva, sol ex')
xlabel('s')
ylabel('V')

%% approssimazione quasi stazionaria

y0_qs = y0(1:2);

[tout_qs yout_qs] = ode15s(@f_es3_1_quasi_staz,t,y0_qs,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,2),'InitialStep',1.e-5),k);

km = (k(2)+k(3))/k(1);
ki = k(5)/k(4);
c1 = k(6)*yout_qs(:,1)./(yout_qs(:,1)+km*(1+yout_qs(:,2)/ki));
c2 = k(6)*yout_qs(:,2)./(yout_qs(:,2)+ki*(1+yout_qs(:,1)/km));

figure(2)
%subplot con s c1 e c2 confrontate tra le rispettive soluzioni
subplot(3,1,1)
plot(tout_qs,yout_qs(:,1),tout_ex,yout(:,1))
title('s(t)')
legend('approx qs','sol ex')

subplot(3,1,2)
plot(tout_qs,c1,tout_ex,yout(:,3))
title('c1(t)')
legend('approx qs','sol ex')

subplot(3,1,3)
plot(tout_qs,c2,tout_ex,yout(:,4))
title('c2(t)')
legend('approx qs','sol ex')

figure(3)
%unico plot con tutte insime con colori diversi e tratteggi
plot(tout_qs,yout_qs(:,1),'--r',tout_ex,yout(:,1),'-r')

hold on
plot(tout_qs,c1,'--g',tout_ex,yout(:,3),'-g')
hold on
plot(tout_qs,c2,'--b',tout_ex,yout(:,4),'-b')
legend('s(t) - approx qs','s(t) - sol ex','c1(t) - approx qs','c1(t) - sol ex','c2(t) - approx qs','c2(t) - sol ex')


V_qs = k(3)*c1;

figure(4)
plot(yout(:,1),V,yout_qs(:,1),V_qs)
axis([0 5 0 inf])
title('confronto V(s) in 2D')
legend('V - sol ex','V - approx qs')
xlabel('s')
ylabel('V')

figure(5)
plot3(yout(:,1),yout(:,2),V,yout_qs(:,1),yout_qs(:,2),V_qs)
legend('V(s,i) competitiva - sol ex','V(s,i) competitiva - approx qs')
title('confronto V(s) in 3D')


%% inibizione allosterica

% t = [0:0.01:100];
y0_al = [5,0,0,0,2];
% k = [102,50,1,26,50,1];

[t yout_al] = ode15s(@f_es3_allo,t,y0_al,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,5),'InitialStep',1.e-5, 'MaxStep', 5),k);

figure(6)
plot(t,yout_al(:,1),'r',t,yout_al(:,2),'m',t,yout_al(:,3),'b',t,yout_al(:,4),'c',t,yout_al(:,5),'g')
legend('s(t)','x(t)','y(t)','z(t)','i(t)')

figure(7)
V_al = k(3)*yout_al(:,2);
plot3(yout_al(:,3),t,V_al)
title('V(s,t)')

y0_al_qs = [5,2];
[t yout_al_qs] = ode15s(@f_es3_allo_qs,t,y0_al_qs,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,2),'InitialStep',1.e-5),k);

figure(8)
plot(t,yout_al_qs(:,1),'--r',t,yout_al(:,1),'r',t,yout_al_qs(:,2),'--b',t,yout_al(:,5),'b')
legend('s(t) - qs','s(t) - intero','i(t) - qs','i(t) - intero')

k1 = k(1);
km1 = k(2);
k2 = k(3);
k3 = k(4);
km3 = k(5);
e0 = k(6);

s = yout_al_qs(:,1);
i = yout_al_qs(:,2);

V_al_qs = (e0*k1*k2*km3*s.*(km1 + km3 + i*k3 + k1*s))./((km3 + i*k3).*(k2*km1 + k2*km3 + km1*km3 + km1^2 + k1^2*s.^2 + i*k3*km1 + k1*k2*s + 2*k1*km1*s + k1*km3*s + i*k1*k3.*s));

% confronto velocità
figure(9)
plot3(yout_al(:,1),yout_al(:,5),V_al, s,i,V_al_qs)
legend('V(s,i) - sistema completo)', 'V(s,i) - approx qs)')
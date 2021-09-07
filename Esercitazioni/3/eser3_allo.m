clear all 
close all
clc

%% inibizione allosterica

t = [0,100];
y0_al = [5,0,0,0,2];
k = [102,50,1,26,50,1];

% sistema completo
[tout_al yout_al] = ode15s(@f_es3_allo,t,y0_al,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,5),'InitialStep',1.e-5, 'MaxStep', 5),k);

figure(1)
plot(tout_al,yout_al(:,1),'r',tout_al,yout_al(:,2),'m',tout_al,yout_al(:,3),'b',tout_al,yout_al(:,4),'c',tout_al,yout_al(:,5),'g')
legend('s(t)','x(t)','y(t)','z(t)','i(t)')
title('5 soluzioni in funzione del tempo, allosterica')
xlabel('t');


figure(2)
V_al = k(3)*yout_al(:,2);
plot3(yout_al(:,3),tout_al,V_al)
title('V(s,t)')
xlabel('s');
ylabel('V');

% approssimazione quasi stazionaria
y0_al_qs = [5,2];
[tout_al_qs yout_al_qs] = ode15s(@f_es3_allo_qs,t,y0_al_qs,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,2),'InitialStep',1.e-5),k);

figure(3)
plot(tout_al_qs,yout_al_qs(:,1),'--r',tout_al,yout_al(:,1),'r',tout_al_qs,yout_al_qs(:,2),'--b',tout_al,yout_al(:,5),'b')
legend('s(t) - qs','s(t) - intero','i(t) - qs','i(t) - intero')
title('confronto s - i allosterica normale e approx qs ')
k1 = k(1);
km1 = k(2);
k2 = k(3);
k3 = k(4);
km3 = k(5);
e0 = k(6);

s = yout_al_qs(:,1);
i = yout_al_qs(:,2);

V_al_qs = (e0*k1*k2*km3*s.*(km1 + km3 + i*k3 + k1*s))./((km3 + i*k3).*(k2*km1 + k2*km3 + km1*km3 + km1^2 + k1^2*s.^2 + i*k3*km1 + k1*k2*s + 2*k1*km1*s + k1*km3*s + i*k1*k3.*s));


% approssimazione all'equilibrio


%% confronto velocità
figure(4)
plot3(yout_al(:,1),yout_al(:,5),V_al, s,i,V_al_qs)
legend('V(s,i) - sistema completo)', 'V(s,i) - approx qs)')
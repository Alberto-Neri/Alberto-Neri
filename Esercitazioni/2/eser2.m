clear all
close all
clc

t = [0,100];
y0 = [5 0 0];

k_1 = [0.5010,500,1,502000,500,2,1]; %coop positiva
k_2 = [102,50,1,26,50,2,1];  %siti indip
k_3 = [1002, 500, 1, 5.02, 500, 2, 1];  %coop negativa

%faccio 3 ode per plottare tutto insieme, se non mi interessa: 
%ne muto 2, cambio i k sopra, muto 2 V
[t1_p yout1] = ode15s(@f_es2_1,t,y0,[],k_1);
[t1_i yout2] = ode15s(@f_es2_1,t,y0,[],k_2);
[t1_n yout3] = ode15s(@f_es2_1,t,y0,[],k_3);

V1 = k_1(3)*yout1(:,2)+k_1(6)*yout1(:,3);
V2 = k_2(3)*yout2(:,2)+k_2(6)*yout2(:,3);
V3 = k_3(3)*yout3(:,2)+k_3(6)*yout3(:,3);

figure(1)
plot(yout1(:,1), V1)
title('Velocità di formazione del composto V(s)')
axis([0, 2, 0, inf])
hold on
plot(yout2(:,1), V2)
hold on
plot(yout3(:,1),V3)
legend('coop pos','siti ind','coop neg');
xlabel('s')
ylabel('V')
%come su quaderno: coop pos sigmoide, indp MM, coop neg piu bassa --> ok

%% approx quasi stazionaria
%su carta devo porre le 2 derivate di c = 0, ricavo c1 e c2 in funzione di
%s, passo tutto alla funzione la cui unica incognita sarà s. Poi uso ode15s
%per risolvere

%DECIDO DI USARE I k DELLA COOP POSITIVA (messi sopra)!!!!!!!!!!
K1_1 = (k_1(2)+k_1(3))/k_1(1);
K2_1 = (k_1(5)+k_1(6))/k_1(4);
K_1 = [k_1, K1_1,K2_1];
s0 = y0(1);


[t_qs s_1] = ode15s(@f_es2_quasi_staz,t,s0,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5),K_1);

V_1 = K_1(7)*s_1.*(K_1(3)*K_1(9)+K_1(6)*s_1)./(K_1(8)*K_1(9)+K_1(9)*s_1+s_1.^2);

figure(2)
plot(s_1,V_1)
hold on
plot(yout1(:,1), V1)
axis([0 2 0 inf])
legend('Approx','Sistema')
title('Velocità di formazione del composto V(s)')
xlabel('s')
ylabel('V')

figure
subplot(3,1,1)
plot(t1_p,yout1(:,1), t_qs,s_1);
title('confronto s(t) - coop positiva');
xlabel('t')
ylabel('s')
legend('Ex','Approx')

e0 = k_1(7);
c1 = K2_1*e0*s_1./(K1_1*K2_1 + K2_1*s_1 + s_1.^2);
c2 = e0*s_1.^2./(K1_1*K2_1 + K2_1*s_1 + s_1.^2);

subplot(3,1,2)
plot(t1_p,yout1(:,2), t_qs,c1);
title('confronto c1(t) - coop positiva');
xlabel('t')
ylabel('c1')
legend('Ex','Approx')

subplot(3,1,3)
plot(t1_p,yout1(:,3), t_qs,c2);
title('confronto c2(t) - coop positiva');
xlabel('t')
ylabel('c2')
legend('Ex','Approx')

%si vede che gli errori diventano accettabili solo dopo una prima fase
%transitoria !!!

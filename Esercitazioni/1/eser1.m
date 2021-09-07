clear all
close all
clc

t = [0,100]'; %con ode15s lascio decidere il passo da lei

% parametri
k = [4.e+6, 25, 15, 1.e-8];
k1 = k(1);
kmeno1 = k(2);
k2 = k(3);
e0 = k(4);
Km = (kmeno1+k2)/k1;

%condizioni iniziali
y0 = [1.e-5,0]; % [s0,c0]
s0 = y0(1);
c0 = y0(2);
[t_ex,yout] = ode15s(@f_es1_1,t,y0,odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13]),k);

figure(1) % riscalo (divido per s0 ed e0) per guardarli sullo stesso grafico
plot(t_ex,yout(:,1)/s0)
title('s(t) e c(t) riscalati');
hold on
plot(t_ex,yout(:,2)/e0)
legend('s(t)', 'c(t)');
xlabel('t');
ylabel('c(t); s(t)');
%vedo che substrato parte da 1 (è riscalato) e poi scende --> ok
%il composto parte da 0, cresce rapidamente, decresce lentamente --> ok
%la costante rapida è quella associata a s, perchè decresce piu veloce

figure(2)
plot(yout(:,1),yout(:,2));
title('spazio delle fasi, soluzione numerica')
xlabel('s(t)');
ylabel('c(t)');
%soluzione numerica, prima cresce, poi scende verticalmente
%come visto su quaderno


%% punto 2, spprossimazione uniforme

su0 = y0(1);
[t_approx,su] = ode15s(@f_es1_2,t,su0,odeset('RelTol',5.e-13 ,'AbsTol',1.e-13),k);

cu = e0*su./(Km+su) -e0*su0/(Km+su0)*exp(-(Km+su0)*k1.*t_approx);

figure
subplot(2,1,1);
plot(t_ex,yout(:,1),t_approx,cu);
legend('esatta','approx');
title('c(t)');

subplot(2,1,2);
plot(t_ex,yout(:,2),t_approx,su);
legend('esatta','approx');
title('s(t)');
%dall'analisi grafica si vede come in fase transitorio l'errore si
%nell'ordine di 10*10^-6 successivamente diminuisce fino a circa 5*10^-6
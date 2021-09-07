%% eser 7
% sistema di HH in 0D

clear all
close all
clc
% 1
figure(1)
T = 100; %faccio fino a 40 che si vede meglio
t = [0,T];

%vettore dei parametri
Iapp = 100;
k = [1 120 36 0.3 115 -12 10.6 Iapp 0]; % Cm, gNa, gK, gL, VNa, Vk, VL, Iapp, num_impulsi
% condizioni iniziali
V0 = 2.7570e-04; % 12, 2.7570e-04, 
m0 = 5.2934e-02;
h0= 5.9611e-01;
n0= 3.1768e-01;
y0 = [V0, m0, h0, n0];
%y0 = [2.7570e-04, 5.2934e-02, 5.9611e-01, 3.1768e-01]; %condizioni resting


[tout,yout] = ode15s(@fHH_es7,t,y0,[],k);
INa = k(2)*yout(:,2).^3.*yout(:,3).*(yout(:,1) - k(5)); %quaderno
IK = k(3)*yout(:,4).^4.*(yout(:,1)-k(6));
IL = k(4)*(yout(:,1)-k(7));
Iion = INa+IK+IL;

gNa = INa./(yout(:,1)-k(5));
gK = IK./(yout(:,1)-k(6));
gL = IL./(yout(:,1)-k(7));

%le tau si calcolano come: 1/(alfa+beta)
%le soluzioni allo stato stazionario: alfa/(alfa+beta)
V = yout(:,1);
alfa_m = 0.1*(25-V).*(exp((25-V)/10)-1).^(-1);
alfa_h = 0.07*exp(-V/20);
alfa_n = 0.01*(10-V).*(exp((10-V)/10)-1).^-1;
beta_m = 4*exp(-V/18);
beta_h = (exp((30-V)/10)+1).^-1;
beta_n = 0.125*exp(-V/80);

tau_n = 1./(alfa_n+beta_n);
tau_m = 1./(alfa_m+beta_m);
tau_h = 1./(alfa_h+beta_h);
n_inf = alfa_n./(alfa_n+beta_n);
m_inf = alfa_m./(alfa_m+beta_m);
h_inf = alfa_h./(alfa_h+beta_h);

figure(1)
subplot(2,2,1);
plot(tout,yout(:,1)); %V sul tempo
title("V(t)");
xlabel('t');
ylabel('V');

subplot(2,2,2);
plot(tout,yout(:,2),'b', tout,yout(:,3),'r',tout,yout(:,4),'g');  %variabili di gating su t
legend("m","h","n");
title("varibili di gating");
xlabel('t');
ylabel('m,h,n');

subplot(2,2,3);
plot(tout,gNa, tout,gK);
legend("gNa","gK");
title("conduttanze");
xlabel('t');
ylabel('gNa,gK');

subplot(2,2,4);
plot(tout,INa,'r', tout,IK,'b',tout,IL,'g', tout,Iion,'k');
legend("INa", "IK", "IL", "Iion");
title("correnti");
xlabel('t');

figure(2)
subplot(2,1,1);
plot(V, tau_n,V, tau_m, V, tau_h)
legend("tau_n","tau_m", 'tau_h');
title("costanti di tempo");
xlabel('V');
ylabel('tau');

subplot(2,1,2);
plot(V, n_inf,V, m_inf, V, h_inf)
legend("n inf","m inf", 'h inf');
title("valori quasi stazionari");
xlabel('V');

%punto 1: Iapp=0, V0=12
%grafici di V, m,n,h e le correnti sembrano realistici --> ok
%grafico costanti di tempo: si vede che tau_m è piu piccola, come visto in
%teoria (0.1/0.2 msec) 
%grafico valori qs: m e h pilotano la corrente di Na e al variare di V una
%cresce e l'altra decresce, n pilota la corrente di K

%punto 5: modifico questa parte cambiando le Iapp e togliendo in fHH il
%controllo sul tempo per l'impulso
% ho verificato le prime e sembra tornino --> ok


% per farlo animato ma non ho voglia di farlo, NON CONTROLLATO: 
% % subplot(2,2,1)
% % al = animatedline;
% % for i=1:length(t)
% %     addpoints(al,t(i),yout(i,1));
% %     drawnow
% % end

% % subplot(2,2,2)
% % al2 = animatedline('LineWidth',1.5,'Color','b');
% % al3 = animatedline('LineWidth',1.3,'Color','r');
% % al4 = animatedline('LineWidth',1,'Color','g');
% % for i=1:length(t)
% %     addpoints(al2,t(i),yout(i,2));
% %     addpoints(al3,t(i),yout(i,3));
% %     addpoints(al4,t(i),yout(i,4));
% %     drawnow
% % end
% 
% % subplot(2,2,4)
% % Na = animatedline('LineWidth',1.5,'Color','r');
% % K = animatedline('LineWidth',1.3,'Color','b');
% % L = animatedline('LineWidth',1,'Color','g');
% % ion = animatedline('LineWidth',1,'Color','k');
% % for i=1:length(t)
% %     addpoints(Na,t(i),INa(i));
% %     addpoints(K,t(i),IK(i));
% %     addpoints(L,t(i),IL(i));
% %     addpoints(ion,t(i),Iion(i));
% %     drawnow
% % end


%% 2
%faccio variare V0 e guardo comportamento del primo grafico (sempre Iapp=0)
%se V>=6.6 parte potenziale

%% 3
%cambio Iapp e vedo che la soglia è tra 65 e 65.5

%% 4 effetto refrattarieta 

Iapp = 100;
DI = [20, 17, 15, 12, 10,8,7, 6, 5,4,3,2,1.5]; %per provarne solo 1 scrivo
%un altro DI = n, poi al limite levo il subplot
%DI = 20;

figure(3)
for di=1:length(DI) 
    k = [1 120 36 0.3 115 -12 10.6 Iapp]; % Cm, gNa, gK, gL, VNa, Vk, VL, Iapp
    % condizioni iniziali
    y0 = [2.7570e-04, 5.2934e-02, 5.9611e-01, 3.1768e-01];
    t_int=0; %tempo dell'intervallo in cui sono
    while (t_int+DI(di) <=T)
        t = [t_int t_int+DI(di)]; %1° ciclo -> [0,20] ..
        k(9) = t_int; %serve per fHH altrimenti controllo sempre a 0.1 ms
        t_int = t_int+ DI(di);  
        [tout,yout] = ode15s(@fHH_es7,t,y0,[],k);
        y0 = [yout(end,1),yout(end,2),yout(end,3),yout(end,4)]; %cambio dati iniziali
        subplot(5,3,di)
        plot(tout,yout(:,1))
        title(DI(di))
        hold on
    end
    
end
%la refrattarietà si sente da un DI di 5, poi con 1.5 c'è un caso
%particolare



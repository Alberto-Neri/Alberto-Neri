%% RMC 0D
clear all
close all
clc

% punto a

T = 100; %100 o 200
t = [0,T];

%vettore dei parametri
Iapp =0; %0; 0.5; 1; 1.5; 2; 2.5; 
b=5;
c=1;
beta=0.1;
epsilon = beta;
delta=1;
gamma=0.25;
e=0.1;
k = [Iapp, b, c, beta, delta, gamma, e]; % Iapp, b, c, beta, delta, gamma, e
k2 = [Iapp, b, c, beta, delta, gamma, epsilon];

% condizioni iniziali
v0 = -0.1;  %0.1 oppure 0.6, cambio io
w0 = -0.1;
y0 = [v0; w0];

[tout,youtFHN] = ode15s(@fhn_es6,t,y0,[],k);
[toutRMC,youtRMC] = ode15s(@f_rmc,t,y0,[],k2);

%nullclines: (potrei usare le yout ma non le vedrei sempre complete)
v = -0.5:0.001:1.5; %da cambiare per vedere qualcosa nei grafici
w = -0.5:0.001:1.5;
nullc1FHN = gamma*w; %w=0
nullc2FHN = 1/c*(b*v.*(v-beta).*(delta-v) + Iapp); %v=0
nullc2_notrasl = 1/c*(b*v.*(v-beta).*(delta-v) ); %v=0 senza Iapp
nullc2RMC = b*(v-beta).*(1-v);
nullc1RMC = gamma*w;

figure()
plot(w,nullc1FHN,'-r',v, nullc2FHN,'-r',w,nullc1RMC,'-b',v,nullc2RMC,'-b');
legend('n1FHN', 'n2FHN', 'n1RMC', 'n2RMC')

%%

% figure()
% subplot(2,1,1)
% plot(tout,yout(:,1))
% xlabel('t')
% ylabel('v')
% title('v(t)')
% 
% subplot(2,1,2)
% plot(tout,yout(:,2));
% xlabel('t')
% ylabel('w')
% title('w(t)')

figure()
subplot(2,1,1)
plot(youtFHN(:,1),youtFHN(:,2),w,nullc1FHN,v, nullc2FHN,v, nullc2_notrasl )
xlabel('v')
ylabel('w')
title('v-w spazio delle fasi FHN')
%legend('andamento','nullc1', 'nullc2');
legend('andamento','nullc1', 'nullc2', 'nullc2 non traslata');

subplot(2,1,2)
plot(youtRMC(:,1),youtRMC(:,2),w,nullc1RMC,v, nullc2RMC)
xlabel('v')
ylabel('w')
title('v-w spazio delle fasi RMC')
%legend('andamento','nullc1', 'nullc2');
%legend('andamento','nullc1', 'nullc2');

    
%analisi punto a) (Iapp = 0)
%vedo che v0 = 0.1 è proprio la soglia, per qualsiasi valore v0<=0.1 
%la traiettoria tende subito al punto di stabilità (0,0). Con valori
%superiori invece si ha che h0(w0)<v0 quindi il potenziale entra in una
%fase di eccitazione che lo porta a tendere al ramo h+, dopo si estingue e
%torna a (0,0)

%analisi punto b) (Iapp != 0)
%con un impulso di corrente la nullc1 si sposta in alto di Iapp, ciò
%comporta uno spostamento del punto di equilibrio (intersezione 2 nullc).
%il punto di equilibrio tornerà ad essere (0,0) passato il tempo di
%stimolazione
%per Iapp = 0.5 --> sia con v0 a 0.1 o 0.6 il potenziale si eccita e tende
%a raggiungere h+ di nullc2 (traslata), dopo che si estingue si riporta su
%nullc2 non traslata e tende a 0,0
%per Iapp = 1 --> come sopra
%lo stesso per le altre, perchè con un tempo di stimolazione di 0.5 il pot non ha
%abbastanza tempo per salire con la nullc2 (traslata). Se infatti provo con
%t=5 in alcuni casi si ha l'effetto che mi sarei aspettato.  




%     figure()
%     for j = 1:length(t)
%         subplot(1,2,1)
%         plot(t(j),yout(j,1),'*r')
%         axis manual
%         axis([0 100 0 1]);
%         drawnow
%         hold on
%         subplot(1,2,2)
%         plot(t(j),yout(j,2),'*b')
%         axis manual
%         axis([0 100 0 1]);
%         drawnow
%         hold on
%     end


%% punto b
%cambio quello che serve sopra
% close all
% 
% T = 200;
% t = 0:0.1:T;
% 
% 
% 
% % condizioni iniziali
% v0 = 0.6;
% w0 = 0;
% y0 = [v0, w0];
% 
% Iapp = 0.5;
% 
% for i = 1:5
%     %vettore dei parametri
%     k = [Iapp*i, 5, 1 0.1, 1, 0.25, 0.1]; % Iapp, b, c, beta, delta, gamma, e
%     
%     [tout,yout] = ode15s(@fhn_es6,t,y0,[],k);
%     figure()
%     subplot(3,1,1)
%     plot(tout,yout(:,1))
%     xlabel('t')
%     ylabel('v')
%     title('v(t)')
%     
%     subplot(3,1,2)
%     plot(tout,yout(:,2));
%     xlabel('t')
%     ylabel('w')
%     title('w(t)')
%     
%     subplot(3,1,3)
%     plot(yout(:,1),yout(:,2))
%     xlabel('v')
%     ylabel('w')
%     title('v-w spazio delle fasi')
%     
% %     figure()
% %     for j = 1:length(t)
% %         subplot(1,2,1)
% %         plot(t(j),yout(j,1),'or')
% %         axis manual
% %         axis([0 200 0 5]);
% %         drawnow
% %         hold on
% %         subplot(1,2,2)
% %         plot(t(j),yout(j,2),'ob')
% %         axis manual
% %         axis([0 200 0 4]);
% %         drawnow
% %         hold on
% %     end
% end

%% HH 0D veloce-lento
% sistema di HH in 0D

clear all
close all
clc
% 1
T = 100; %faccio fino a 40 che si vede meglio
t = [0,T];

%vettore dei parametri
Iapp = 100;
k = [1 120 36 0.3 115 -12 10.6 Iapp 0]; % Cm, gNa, gK, gL, VNa, Vk, VL, Iapp, num_impulsi
Cm=k(1);
gNa =k(2); %sono massimali
gK=k(3);
gL=k(4);
VNa = k(5);
VK=k(6);
VL=k(7);
Iapp = k(8);
t0 = k(9);

% condizioni iniziali
V0 = 0.1; % 12, 2.7570e-04, 
m0 = 5.2934e-02;
h0= 5.9611e-01;
n0= 0.3; %3.1768e-01
y0 = [V0, m0, h0, n0];
%y0 = [V0, 5.2934e-02, 5.9611e-01, 3.1768e-01]; %condizioni resting


[tout,yout] = ode15s(@fHH_es7,t,y0,[],k); %simulo tutto anche m,h !!
Vout = yout(:,1);
nout = yout(:,4);
% INa = k(2)*yout(:,2).^3.*yout(:,3).*(yout(:,1) - k(5)); %quaderno
% IK = k(3)*yout(:,4).^4.*(yout(:,1)-k(6));
% IL = k(4)*(yout(:,1)-k(7));
% Iion = INa+IK+IL;
% 
% gNa = INa./(yout(:,1)-k(5));
% gK = IK./(yout(:,1)-k(6));
% gL = IL./(yout(:,1)-k(7));

%le tau si calcolano come: 1/(alfa+beta)
%le soluzioni allo stato stazionario: alfa/(alfa+beta)
% V = yout(:,1);
[V,n] = meshgrid(-20:0.1:120, 0:0.01:1); %definisco i punti su cui voglio vedere le nullc  
%se non faccio la grid avrei 2 vettori di dimensioni diverse che poi si
%moltiplicano in nullc_V perciò non verrà mai!!!

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

nullc_V = gNa*m_inf.^3.*(0.8-n).*(V-VNa) + gK*n.^4.*(V-VK) + gL*(V-VL); %scrivo solo termine noto!!
nullc_n = n_inf;

figure(1)
subplot(2,1,1)
contour(V,n,nullc_V, [0 0])
hold on
plot(V(1,:),nullc_n(1,:))
xlabel('V')
ylabel('n')
plot(Vout, nout);
legend('nullc-V','nullc-n','traiett')

subplot(2,1,2)
plot(tout, Vout);
xlabel('t')
ylabel('V')
%P è a t=4.207 (quando v raggiunge il minimo)

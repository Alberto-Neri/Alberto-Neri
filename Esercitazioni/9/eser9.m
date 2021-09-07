%% Modello HH 1D: propagazione potenziale d'azione lungo fibra

close all
clear all
clc

T = 20;  %serve tempo >1 se no non si vede niente
k=0.010; %sarebbe il tau --> variazione nel tempo
tempo=0:k:T;
nt=T/k;

%vettore dei parametri
Iapp = 100;
Cm=1;
gNa =120; %sono massimali
gK=36;
gL=0.3;
VNa = 115;
VK=-12;
VL=10.6;
sigma = 0.001; %lo prendo dalla 8 (FHN 1D)

% condizioni iniziali di resting
v0 = 2.7570e-04;  
m0 = 5.2934e-02;
h0= 5.9611e-01;
n0= 3.1768e-01;
%y0 = [2.7570e-04, 5.2934e-02, 5.9611e-01, 3.1768e-01]; %condizioni resting


for nn=100

  t=0;
  tindex=1;
  hh=1/(nn+1); %variazione nello spazio
  x=[0:hh:1]';
  
  %condizioni iniziali v0, m0, h0, n0 di resting
  v(:,tindex)=v0.*ones(length(x),1); 
  m(:,tindex)=m0*ones(length(x),1);
  h(:,tindex)=h0*ones(length(x),1);
  n(:,tindex)=n0*ones(length(x),1);
  
  A=-2*diag(ones(nn,1))+diag(ones(nn-1,1),1)+diag(ones(nn-1,1),-1);
  A = [-2 2 zeros(1,nn); [1; zeros(nn-1,1)], A, [zeros(nn-1,1); 1]; zeros(1,nn) 2 -2]; %per condizioni di neumann --> ghost point
  
  AA = Cm*eye(nn+2)/k - sigma/hh^2*A; %per conti quaderno
 
  bb=hh^2*zeros(nn+2,1);
  Iapp = zeros(length(x),length(tempo));
  
  %Iapp(find(x<=0.04)) = 10
  %trovo modo per definire Iapp=10 quando x<=0.04 && t<1 
  %vale per estremo dx
  for i=1:length(x)
      for j=1:length(tempo)
          if (x(i) <= 0.04 && tempo(j)<=1) 
              Iapp(i,j) = 100;
          end
%           if (x(i) >= x(end)-0.04 && tempo(j)<=1) %impulso da sx
%               Iapp(i,j) = 100;
%           end
          
      end
  end
  
  for kk=1:nt
          
    v_j = v(:,tindex); %per migliore lettura, sarebbe v_t
    m_t = m(:,tindex);
    h_t = h(:,tindex);
    n_t = n(:,tindex);
    
    bb= Cm*v_j/k -(gNa*m_t.^3.*h_t.*(v_j-VNa) + gK*n_t.^4.*(v_j-VK) + gL*(v_j-VL)) + Iapp(:,tindex);
    
    vt=AA\bb; %sto calcolando v_t+1
    
    %calcolo alfa, beta per il t corrente, quindi uso v_j
    am = 0.1*(25-v_j).*(exp((25-v_j)/10)-1).^(-1);
    ah = 0.07*exp(-v_j/20);
    an = 0.01*(10-v_j).*(exp((10-v_j)/10)-1).^-1;
    bm = 4*exp(-v_j/18);
    bh = (exp((30-v_j)/10)+1).^-1;
    bn = 0.125*exp(-v_j/80);
    
    %calcolo m_t+1
    m(:,tindex+1) = m_t + k*(am.*(1-m_t)-bm.*m_t); %conti quaderno   
    %calcolo h_t+1
    h(:,tindex+1) = h_t + k*(ah.*(1-h_t)-bh.*h_t); %conti quaderno
    %calcolo n_t+1
    n(:,tindex+1) = n_t + k*(an.*(1-n_t)-bn.*n_t); %conti quaderno
    
    tindex=tindex+1;
    v(:,tindex)=vt; %aggiungo vt come vettore colonna a v_j al ciclo dopo, e diventerà matrice
    t = t+k;
    
    %%%%%%MEGLIO FAR ANDARE 1 FIGURA ALLA VOLTA%%%%%%%%%%%
    figure(1)
    plot(x,vt,'-o');  
    title(num2str(t))
    axis manual
    axis([0 1 -20 120]); 
    xlabel('x');
    ylabel('v(x,t)');
    drawnow
    
%     figure(2)
%     subplot(2,1,1)
%     plot(x,m_t,x,h_t,x,n_t);
%     title(num2str(t))
%     xlabel('x');
%     ylabel('var gating');
%     legend('m','h','n');
%     
%     %calcolo le correnti di Na K e ionica
%     INa = gNa*m_t.^3.*h_t.*(vt - VNa); %quaderno
%     IK = gK*n_t.^4.*(vt-VK);
%     IL = gL*(vt-VL);
%     Iion = INa+IK+IL;
%     subplot(2,1,2)
%     plot(x,INa,x,IK,x,Iion);
%     title(num2str(t))
%     xlabel('x');
%     ylabel('correnti');
%     legend('INa','IK','Iion');
%     drawnow

  end

end
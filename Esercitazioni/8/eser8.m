%% Modello FHN 1D: propagazione potenziale d'azione lungo fibra

close all
clear all
clc

T = 20;  %serve tempo >1 se no non si vede niente
k=0.010; %sarebbe il tau --> variazione nel tempo
tempo=0:k:T;
nt=T/k;
sigma = 0.001;
b = 5;
c=1;
beta = 0.1;
delta = 1;  %Ipotesi nagumo
gamma = 0.25;
e = 0.1;
% index=0;
%for n=[40 80 160 320];

for n=100

  t=0;
  tindex=1;
  h=1/(n+1); %variazione nello spazio
  x=[0:h:1]';
  
  %condizioni iniziali v0, w0  ??? per ora lascio a zero
  v(:,tindex)=0.*ones(length(x),1); 
  w(:,tindex)=0.*ones(length(x),1);
  
  A=-2*diag(ones(n,1))+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
  A = [-2 2 zeros(1,n); [1; zeros(n-1,1)], A, [zeros(n-1,1); 1]; zeros(1,n) 2 -2]; %per condizioni di neumann --> ghost point
  
  AA = eye(n+2)/k - sigma/h^2*A; %per conti quaderno
 
  bb=h^2*zeros(n+2,1);
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
          
    v_j = v(:,tindex); %per migliore lettura, sarebbe vt
    
    bb= v_j/k +( b*v_j.*(v_j-beta).*(delta-v_j) )- c*w(:,tindex) +Iapp(:,tindex);
    
    vt=AA\bb; %sto calcolando v_t+1
    
    %calcolo w_t+1
    w(:,tindex+1) = w(:,tindex)*(1-k*e*gamma) + k*e*v_j; %conti quaderno  
    
    figure(1)
    tindex=tindex+1;
    v(:,tindex)=vt; %aggiungo vt come vettore colonna a v_j al ciclo dopo, e diventerà matrice
    t = t+k;
     
    plot(x,vt,'-o');
    title(num2str(t))
    axis manual
    axis([0 1 -1 3.5]); 
    xlabel('x');
    ylabel('v(x,t)');
    drawnow
  end
%   figure(2)
%   plot(tt,err);
%   hold on
%   xlabel('t');
%   ylabel('err');
end
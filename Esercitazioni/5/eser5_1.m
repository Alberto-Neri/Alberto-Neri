%% esercitazione 5.1 kolmogorov
clear all;
close all;
clc

T = 10;  %serve tempo >1 se no non si vede niente
k=0.010;
tempo=0:k:T;
nt=T/k;
sigma = 0.001;
b = 5;
beta = 0.1;
delta = 1;  %Ipotesi nagumo
% index=0;
%for n=[40 80 160 320];

for n=[40 80];
  clear u x A B bb tt err uex
  t=0;
  tindex=1;
  h=1/(n+1);
  x=[0:h:1]';
  u(:,tindex)=zeros(length(x),1); %condizione iniziale
  A=-2*diag(ones(n,1))+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
  A = [-2 2 zeros(1,n); [1; zeros(n-1,1)], A, [zeros(n-1,1); 1]; zeros(1,n) 2 -2]; %per condizioni di neumann --> ghost point
  
  AA = eye(n+2)/k - sigma/h^2*A; %per conti quaderno
 
  bb=h^2*zeros(n+2,1);
  Iapp = zeros(length(x),length(tempo));
  
  %Iapp(find(x<=0.04)) = 10
  %trovo modo per definire Iapp=10 quando x<=0.04 && t<1 
  for i=1:length(x)
      for j=1:length(tempo)
          if (x(i) <= 0.04 && tempo(j)<=1) 
              Iapp(i,j) = 10;
          end
      end
  end
  
  for kk=1:nt
    
    %la prima è fisher, la seconda Nagumo, muto quella che non mi serve
    u_j = u(:,tindex); %u_t attuale
    %%bb = u_j/k + b*u_j.*(1-u_j)+Iapp(:,tindex);
    bb= u_j/k +( b*u_j.*(u_j-beta).*(delta-u_j) )+Iapp(:,tindex);
    
    ut=AA\bb; %sto calcolando u_t+1
    
    figure(1)
    tindex=tindex+1;
    u(:,tindex)=ut; %aggiungo ut come vettore colonna a u_j al ciclo dopo, e diventerà matrice
    t = t+k;
  
    plot(x,ut,'-o');
    title(num2str(t))
    axis manual
    %axis([0 1 0 5]);  %per kolmogorov
    axis([0 1 0 2]);  %per nagumo.. non bello
    xlabel('x');
    ylabel('u(x,t)');
    drawnow
  end;
%   figure(2)
%   plot(tt,err);
%   hold on
%   xlabel('t');
%   ylabel('err');
end;

%figure(3)
% semilogy(err,'-o');
% mesh(tt,x,u)
% xlabel('t');
% ylabel('x');
% zlabel('u');


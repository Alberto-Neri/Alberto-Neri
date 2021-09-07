clear all;
close all;
clc

T = 5;  %serve tempo >1 se no non si vede niente
k=0.010;
nt=T/k;
sigma = 0.001;
b = 5;
beta = 0.1;
delta = 1;
n=80;
tindex=1;
h=1/(n+1);
x=[0:h:1]';
t=[0:k:T]';

uu =0;
u(:,tindex) = zeros(length(x),1);

A=2*diag(ones(n+2,1))-diag(ones(n+1,1),1)-diag(ones(n+1,1),-1);
A(1,2)=-2;
A(end, end-1) = -2;

A = A*(1/h^2);
AA = eye(n+2)/k + sigma*A;
bb=h^2*zeros(n+2,1);
Iapp = zeros(length(x),length(t));

for i=1:length(x)
      for j=1:length(t)
          if (x(i) <= 0.04 && t(j)<=1) 
              Iapp(i,j) = 2;
          end
      end
end
  
for kk=1:nt
    
    
   u_j = u(:,tindex);
   t_noto =  b*u_j.*(u_j-beta).*(delta-u_j);
   bb= u_j/k +(t_noto )+Iapp(:,tindex);
   tindex = tindex+1;
   
   u(:,tindex) = AA\bb;
   
   figure(1)
    
   
  
    plot(x,u(:,tindex),'-o');
    title(num2str(t(kk)))
    axis manual
    %axis([0 1 0 5]);  %per kolmogorov
    axis([0 1 0 2]);  %per nagumo.. non bello
    xlabel('x');
    ylabel('u(x,t)');
    drawnow



end
figure(2)
  plot(t,u(40,:));
  
  xlabel('t');
  ylabel('ut');

  
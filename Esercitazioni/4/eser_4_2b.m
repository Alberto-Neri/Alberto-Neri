%% 4.2 B) con condizioni di Neumann
clear all;
close all;

T = 0.4;
k=0.001;
nt=T/k;
% index=0;
%for n=[40 80 160 320];
for n=[40 80];
  clear u x A B b tt err uex
  t=0;
  tindex=1;
  h=1/(n+1);
  x=[0:h:1]';
  u(:,tindex)=100*sin(pi*x);
  A=2*diag(ones(n+2,1))-diag(ones(n+2-1,1),1)-diag(ones(n+2-1,1),-1); %cambio A
  A(1,2)=-2;
  A(end,end-1)=-2;
  B = eye(n+2)+(k/h^2)*A; %deve essere della stessa dimensione di x
  b=zeros(n+2,1);
  for kk=1:nt
    b = u(:,tindex);
    ut=B\b;
    %ut=[0; ut; 0]; %condizione di dirichlet che tolgo
    figure(1)
    tindex=tindex+1;
    u(:,tindex)=ut;
    t = t+k;
    uex=100*exp(-pi^2*t).*sin(pi*x); %non è piu la soluzione esatta
%    err(tindex) = norm(ut-uex,inf);
    tt(tindex) = t;
    plot(x,ut,'-o',x,uex,'k--');
    %plot(x,ut);
    %legend('ut','uex');
    title(num2str(t))
    axis manual
    axis([0 1 0 100]);
    xlabel('x');
    ylabel('u(x,t)');
    drawnow
  end;

%   plot(tt,err);
%   hold on
%   xlabel('t');
%   ylabel('err');
end;

% figure(3)
% semilogy(err,'-o');
% mesh(tt,x,u)
% xlabel('t');
% ylabel('x');
% zlabel('u');


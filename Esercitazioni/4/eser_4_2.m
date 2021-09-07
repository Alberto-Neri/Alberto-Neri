% eser4

clear all
close all

%% 4.2 a) Discretizzare l'equazione del calore in 1D in spazio
%con dati al bordo Dirichlet

%u_t = u_xx

X=100;
h= 1/(X+1); %passo di discretizzazione (SPAZIO)
x=0:h:1; %griglia (SPAZIO)

T = 1;
k = 0.01; %passo di discretizzazione (TEMPO)
t=0:k:T;  %griglia (TEMPO)

%costruisco matrice A
n = length(x);
A= zeros(n);
for i=1:n 
    for j=1:n
        if (i==j)
            A(i,j)=2;
            if (i==1)
                A(i,j+1)=-1;
            elseif (i==n)
                A(i,j-1)=-1;
            else 
               A(i,j+1)=-1;
               A(i,j-1)=-1;
        
            end 
        end 
       
    end 
end

%matrice identita
I = eye(n);

u = zeros(length(x),length(t));  %righe: spazio, colonne: tempo
u(1,:) = 0; %boundary condition
u(:,1) = 100*sin(3.14*x);
u(X,:) = 0;

u_ex = zeros(n,1);
err = zeros(length(t)-1,1);
%ora devo trovare gli altri valori con eulero implicito
%avendo definito u come matrice ho bisogno di 2 indici, i:righe(spazio) e
%j:colonne(tempo)

for j=1:length(t)-1  
     %su ciascuna riga faccio eulero implicito per trovare u di j+1
     u(:,j+1) = (I + k/h^2*A)\u(:,j); %definisco la colonna j+1
     for i=1:length(x)
        u_ex(i) = 100*exp(-3.14^2*t(j))*sin(3.14*x(i)); % verifico soluzione esatta
     end
     err(i) = norm(u(:,j+1) - u_ex, inf); %calcolo l'errore su vettore
end %dovrei aver riempito u

%grafico
surf(t,x,u);

figure(2)
plot(x,err);
title('errore'); %è tutto zero --> ho verificato la soluzione esatta


%% b) con dati al bordo Neumann ???



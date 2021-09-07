%Esercitazione 4
%Discretizzazione dell'equazione di Poisson in 1D

n=100;
h= 1/(n+1); %passo di discretizzazione
x=0:h:1; %griglia

u_esatta= -x.^2/2 -3*x/2 +1;
u= zeros(n,1);
u(1)=1;
u(n)=-1;
senx = sin(x);

%costruisco matrice b
b= zeros(n,1);
for i= 1:n
    if (i==1)
        b(i)= senx(i) + (1/h.^2)*u(1);
    else if  (i==n)
        b(n)= senx(i) + (1/h.^2)*u(n);
    else b(i)=1;
    
        end 
    end 
end 

%costruisco matrice A
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

A=(1/h.^2)*A;
u= A\b; %soluzione sist dinamico 
u_approx=[1;u;-1]; %incollo le condizioni al bordo 

figure(1)
plot(x,u_approx,'*r');
hold on 
plot(x,u_esatta,'b');
title('Confronto soluzione esatta e approx');
legend('approx','esatta')

errore= u_esatta - u_approx';
Einf = norm(u_esatta -u_approx', Inf);

figure(2)
loglog(x,errore);
title('Errore')
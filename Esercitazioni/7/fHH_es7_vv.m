function dy = fHH_es7_vv(t,y,k)
    dy = zeros(4,1); % V m h n
    dV = dy(1);
    dm = dy(2);
    dh = dy(3);
    dn = dy(4);
    
    Cm=k(1);
    gNa =k(2); %sono massimali
    gK=k(3);
    gL=k(4);
    VNa = k(5);
    VK=k(6);
    VL=k(7);
    Iapp = k(8);
    t0 = k(9);
    
    %da mutare tutto if per pt.5 (stimolazione continua)
    if t>0.1+t0 %+t0 serve per pt.4 altrimenti quando lancio la ode su un
                %intervallo diverso da [0 ..] non avrò mai l'impulso di
                %corrente visto che controllerei solo su 0.1
        Iapp =0;
            
    end
    
%     if (t<2.2 || t>2.3)
%         Iapp =0;        
%     end
    
    V = y(1);
    m = y(2);
    h = y(3);
    n = y(4);
    
    am = 0.1*(25-V)*(exp((25-V)/10)-1)^(-1);
    ah = 0.07*exp(-V/20);
    an = 0.01*(10-V)*(exp((10-V)/10)-1)^-1;
    bm = 4*exp(-V/18);
    bh = (exp((30-V)/10)+1)^-1;
    bn = 0.125*exp(-V/80);
    

dy(1) = 1/Cm*(-(gNa*m^3*h*(V-VNa) + gK*n^4*(V-VK) + gL*(V-VL)) + Iapp);
dy(2) = am*(1-m) - bm*m;
% dy(3) = h;  %voglio che restino costanti!!!
% dy(4) = n;

end
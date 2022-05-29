Pt_dBm=2; %Input transmitted power in dBm
Gt_dBi=3; %Gain of the Transmitted antenna in dBi
Gr_dBi=3; %Gain of the Receiver antenna in dBi
f=2.9e9; %Transmitted signal frequency in Hertz
d0=1; %assume reference distance = 1m
d=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]; %Array of distances to simulate
Pr_test=[-74.9 -69.4 -70.6 -71.9 -73.6 -76.2 -81.6 -88.2 -75.8 -93.5 -79.2 -83.1 -82.7 -83.2 -82.5 -80.2 -90.1 -85.7 -83.6 -84.9 -83.8 -85.5 -84.6];
L=1; %Other System Losses, No Loss case L=1
sigma=2;%Standard deviation of log Normal distribution (in dB)
n=2; % path loss exponent

Pr_shadow = -1 * Generate_PL_shadow(d,d0);

Pr_Friss = Generate_Pr_friis(Pt_dBm,Gt_dBi,Gr_dBi,Generate_PLdo(d),L);

display(Pr_shadow)
display(Pr_Friss)
display(Pr_test)

PL = logNormalShadowing(Pt_dBm,Gt_dBi,Gr_dBi,f,d0,d,L,sigma,n)

display(PL)

figure;plot(d,PL,'b');hold on;
plot(d,Pr_Friss,'r');grid on;
plot(d,Pr_test,'g');grid on;
xlabel('Distance (m)'); ylabel('P_r (dBm)');
title('Log Normal Shadowing Model');
legend('Log normal shadowing','Friss model', 'Recibido');

function [PL,Pr] = logNormalShadowing(Pt,Gt,Gr,f,d0,d,L,sigma,n)
    lambda = (3*10^8)/f;
    K = 20*log10(lambda/(4*pi)) - 10*n*log10(d0) - 10*log10(L);
    X = sigma*randn(1,numel(d)); 

    PL = Gt + Gr + K - 10*n*log10(d/d0) - X ;
    Pr = Pt + PL;
end

function Pr_shadow = Generate_Pr_friis(Pt,Gt,Gr,Pl,L)
    size = length(Pl);
    for i=1:size
        Pr_shadow(i) = Pt + Gt + Gr - L - Pl(i);
    end
    
end

function PL_shadow = Generate_PL_shadow(array,d0)
    c = 3e8;
    f1 = 2.97e9;
    lambda = c/f1;
    po = 20*log((4*pi*d0/lambda));
    X = 2*randn(1,numel(array));
    size = length(array);
    PL_shadow = po + (20*log10(array/d0)) + X;
end

function PLdo = Generate_PLdo(array)
    size = length(array);
    c = 3e8;
    f1 = 2.97e9;
    lambda = c/f1;
    for i=1:size
        PLdo(i) = 10*log10((4*pi*array(i)/lambda)^2);
    end
end

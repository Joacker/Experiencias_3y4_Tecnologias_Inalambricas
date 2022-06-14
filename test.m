Pt_dBm=0; %Input transmitted power in dBm
Gt_dBi=3; %Gain of the Transmitted antenna in dBi
Gr_dBi=3; %Gain of the Receiver antenna in dBi
f=2.9e9; %Transmitted signal frequency in Hertz
d0=1; %assume reference distance = 1m
d=[1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 9 9.5 10 11 12 13 14 15 16 17 18 19 20 21]; %Array of distances to simulate
Pr_test=[-74.9; -42; -69.4; -49; -70.6; -53.5; -71.9; -50; -73.6; -53.5; -76.2; -50; -81.6; -51.5; -88.2; -75.8; -93.5; -79.2; -83.1; -82.7; -83.2; -82.5; -80.2; -90.1; -85.7; -83.6; -84.9; -83.8; -85.5];
L=1; %Other System Losses, No Loss case L=1
sigma=2;%Standard deviation of log Normal distribution (in dB)
n=2; % path loss exponent
X = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 9 9.5 10 11 12 13 14 15 16 17 18 19 20 21]
Y = [-82.4; -85.1; -90.4; -90.9; -90.5; -92.8; -95.6; -94.5; -91.6; -91.8; -93.2; -96.6; -97.5; -101.3; -105.8; -101.9; -100.3; -100.1; -98.6; -99.1; -100.6; -95.5; -98; -100; -102.3; -100.2; -102.7; -102.7; -110.8]

Pr_shadow = Generate_PL_shadow(d,d0);

Pr_Friss = Generate_Pr_friis(Pt_dBm,Gt_dBi,Gr_dBi,Generate_PLdo(d),L);

%display(Pr_shadow)
%display(Pr_Friss)
%display(Pr_test)
u = 21:1:49
newD = [d,u]

%Extracting the midponts
mid_Y1 = ExtractMidPoints(Y);
mid_Pr_test = ExtractMidPoints(Pr_test);

% EXP3
exp3 = [mid_Y1(:) Y]' ;
composeExp3 = exp3(:);

% EXP4
exp4 = [mid_Pr_test(:) Y]' ;
composeExp4 = exp4(:);

% EXP3 + EXP4
newArray = [Pr_test Y]' ;
c = newArray(:);
c2 = [c]';

PL = logNormalShadowing(Pt_dBm,Gt_dBi,Gr_dBi,f,d0,d,L,sigma,n);

% Test Mid
Mid = MidPoints(Y,Pr_test);

% Pasar a veces los dbm
v1 = ParseVeces(Pr_test);
v2 = ParseVeces(Y);
v3 = ParseVeces(c2);

% Pasar a db las veces
db1 = ParsedB(v1);
db2 = ParsedB(v2);
db3 = ParsedB(v3);

% muestras con ubicaciones
M1 = [d(:),v1(:)];
M2 = [d(:),v2(:)];
M3 = [newD(:), v3(:)];

% Muestras concatenadas y generacion de rectas
y1 = MeanSquareError(M1);
y2 = MeanSquareError(M2);
y3 = MeanSquareError(M3);

% Estadisticas asociadas
err1 = y1 - M1(:,2);
err2 = y2 - M2(:,2);
err3 = y3 - M3(:,2);

test = [newD(:), err3(:)];

% Paso a db en error
nerr1 = ParsedB(err1);
nerr2 = ParsedB(err2);
nerr3 = ParsedB(err3);

% Último ajuste de métodos cuadrados 
a1 = MetodosCuadrados(d,ParsedB(y1));
a2 = MetodosCuadrados(d,ParsedB(y2));
a3 = MetodosCuadrados(newD,ParsedB(y3));

% Último ajuste para pasar los imaginarios a reales
werr1 = real(nerr1);
werr2 = real(nerr2);
werr3 = real(nerr3);

% Media y desviación estandar exp3
media1 = mean(werr1);
desv1 = std(werr1);

% Media y desviación estandar exp4
media2 = mean(werr2);
desv2 = std(werr2);

% Media y desviación estandar exp3 + exp4
media3 = mean(werr3);
desv3 = std(werr3);

% Medias
display(media1)
display(media2)
display(media3)

% Desviaciones Estandar
display(desv1)
display(desv2)
display(desv3)

% Plots
figure;
cdfplot(werr1(:));hold on;
cdfplot(werr2(:));grid on;
cdfplot(werr3(:));grid on;
xlabel('dB'); ylabel('Probabilidad');
legend('Exp3','Exp4','Exp3 + Exp4');

figure;
cdfplot(werr1(:));hold on;
xlabel('dB'); ylabel('Probabilidad');
legend('Exp3');

figure;
cdfplot(werr2(:));hold on;
xlabel('dB'); ylabel('Probabilidad');
legend('Exp4');

figure;
cdfplot(werr3(:));hold on;
xlabel('dB'); ylabel('Probabilidad');
legend('Exp3 + Exp4');

figure;
plot(d,a1,'b');hold on;
plot(d,nerr1,'o');grid on;
xlabel('Distance (m)'); ylabel('dB');
title('MeanSquareError');
legend('Exp3','error exp3');

figure;
plot(d,a2,'b');hold on;
plot(d,nerr2,'o');grid on;
xlabel('Distance (m)'); ylabel('dB');
title('MeanSquareError');
legend('Exp4','error exp4');

figure;
plot(newD,a3,'b');hold on;
plot(newD,nerr3,'o');grid on;
xlabel('Distance (m)'); ylabel('dB');
title('MeanSquareError');
legend('Exp3 + Exp4','error exp3 + exp4');

figure;
%plot(d,Pr_test,'g');hold on;
%plot(d,MetodosCuadrados(d,Pr_test),'c');grid on;
plot(newD,c2,'k');hold on;
plot(newD,MetodosCuadrados(newD,c2),'c');grid on;
%plot(d,Y,'c');grid on;
%plot(d,MetodosCuadrados(d,Y),'k');grid on;
xlabel('Distance (m)'); ylabel('Pr (dBm)');
title('Aplicación de Metodos de Cuadrados');
legend('Exp 3 + Exp 4','Metodos Cuadrados');

function [y] = ParsedB(x)
    size = length(x);
    for i=1:size
        y(i) = 10*log10(x(i));
    end
end

function [y] = ParseVeces(x)
    size = length(x);
    for i=1:size
        y(i) = 10^((x(i))/10);
    end
end

function [y] = MeanSquareError(M)
    n = length(M(:,1));

    x_sum = sum(M(:,1));
    y_sum = sum(M(:,2));

    x2_sum = sum(M(:,1).^2);
    xy_sum = sum(M(:,1).*M(:,2));

    a = (n*xy_sum - (x_sum*y_sum))/(n*x2_sum - (x_sum^2));
    b = (y_sum/n) - (a.*(x_sum/n));

    y = M(:,1).*a + b;
end

function [MidArray] = ExtractMidPoints(x)
    size = length(x);
    for i=1:size
        if(i <= size-1)
           MidArray(i) = x(i) + x(i+1) / 2;
        else
             MidArray(i) = x(i);
        end
    end
end

function [MidArray] = MidPoints(x,y)
    size = length(x);
    for i=1:size
        MidArray(i) = (y(i) + x(i)) / 2;
    end
end

% y = mx + b
% metodos de mínimos cuadrados
function [objectiveFunction] = MetodosCuadrados(x,y)
    size = length(x);
    xy = zeros(size,1);
    sumX = 0;
    sumY = 0;
    sumXY = 0;
    sumxpow2 = 0;
    xpow2 = zeros(size,1);
    for i=1:size
        xy(i) = x(i)*y(i);
        xpow2(i) = x(i)*x(i);
        sumX = sumX + x(i);
        sumY = sumY + y(i);
        sumXY = sumXY + xy(i);
        sumxpow2 = sumxpow2 + xpow2(i);
    end
    m = (sumXY - ((sumX * sumY)/size)) / (sumxpow2 - ((sumX)^2)/size);
    b = (sumY / size) - m*(sumX / size);
    display(m)
    for i=1:size
        objectiveFunction(i) = m*x(i) + b;
    end
end

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
    PL_shadow = po + (20*log10(array/d0));
end

function PLdo = Generate_PLdo(array)
    size = length(array);
    c = 3e8;
    f1 = 2.97e9;
    lambda = c/f1;
    for i=1:size
        PLdo(i) = 20*log10((4*pi*array(i)/lambda));
    end
end


function Yprima = Intercalate(Exp1,Exp2)
    size = length(Exp1);
    for i=1:size
        if(mod(i,2) == 0)
           Yprima(i) = Exp1(i);
        else
           Yprima(i) = Exp2(i);
        end
    end
end
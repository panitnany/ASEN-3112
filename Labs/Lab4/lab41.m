% --- 
% 
% Lab 4 Structures
% 
% ---
%% housekeepin
clear; clc; close;
%% disp in inches // voltages in V
dispSquare1 = [ 1/8 , 1/4 , 3/8 , 1/2 , 5/8 , 6/8 ];
voltSquare1 = [ .555  .586 .532 .488 .417 .375 ];

dispRect1 = [1/16 2/16 3/16 4/16 5/16 6/16 7/16 8/16 9/16 10/16 11/16 ...
    12/16 13/16 14/16 15/16 16/16 17/16 18/16];
voltRect1 = [.323 .334 .338 .338 .339 .339 .339 .339 .339 ... 
    .339 .338 .336 .333 .330 .330 .327 .319 .310];

% julians data
dispSquare2 = [ 0 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00];
voltSquare2 = [.268 .575 .5821 .563 .528 .445 .379 .315 .284];

dispRect2 = [ 0 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00];
voltRect2 = [ .3185 .3523 .3523 .3528 .3523 .3465 .331 .310 .295];

% cav data 

%% converstion to lbs from V
conv = 1000/2.37;
voltSquare1 = voltSquare1*conv;
voltRect1 = voltRect1*conv;
voltSquare2 = voltSquare2*conv;
voltRect2 = voltRect2*conv;

%% ---- BREAK POINT IN CODEE ----


L = 10.5; %inches
L = L + 0.4; %connector pieces

E_square = 10000000;
E_rect = 10000000;

I_square = .0003051757813;
I_rect = .0001627604167;

stress_square = 35000;
stress_rect = 35000;

strain_square = stress_square/E_square;
strain_rect = stress_rect/E_rect;

P_cr_square = (pi^2*E_square*I_square)/L^2;
P_cr_rect = (pi^2*E_rect*I_rect)/L^2;


delt = linspace(0,2,10000);
x = L/2;
v = @(del) del*sin((pi*x)/L);
k_square = @(del) del*(pi^2/L^2)*sin((pi*x)/L)*0.125;
k_rect = @(del) del*(pi^2/L^2)*sin((pi*x)/L)*0.0625;

for i = 1:length(delt)
    vx(i) = v(delt(i));
    k_s_x = k_square(delt);
    k_r_x = k_rect(delt);
end
tol = 0.00001;
buck_s = delt(find(abs(k_s_x - double(strain_square))<tol,1));
buck_r = delt(find(abs(k_r_x - double(strain_rect))<tol,1));

%% analytical prediction
delta = linspace(0,2,200);
Ps = P_cr_square * ( 1 + (pi^2 /(8*L^2))*delta.^2);
Pr = P_cr_rect * ( 1 + (pi^2 /(8*L^2))*delta.^2);

%% Ploots
figure
plot(dispSquare1,voltSquare1,'o','LineWidth',1.2)
hold on
df = linspace(100,300,100);
yeet = zeros(1,length(df)) + buck_s;
plot(delta,Ps,'LineWidth',1.3)
plot(yeet,df,'--g','LineWidth',1.3)
title('Load vs. Horizontal Deflection, Hollow Square Tube')
xlim([0,2])
ylim([100,300])
% crit line
lb = zeros(1,length(0:.25:2)) + P_cr_square;
plot(0:.25:2,lb,'--k','LineWidth',1.3)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','E')

figure
plot(dispRect1,voltRect1,'o','LineWidth',1.2)
title('Load vs. Horizontal Deflection, Rectangular Bar')
xlim([0,2])
ylim([120,160])
hold on
df = linspace(120,160,100);
yeet = zeros(1,length(df)) + buck_r;
plot(delta,Pr,'LineWidth',1.3)
plot(yeet,df,'--g','LineWidth',1.3)
lb = zeros(1,length(0:.25:2)) + P_cr_rect;
plot(0:.25:2,lb,'--k','LineWidth',1.3)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force')


figure
plot(dispSquare2,voltSquare2,'o','LineWidth',1.2)
title('Load vs. Horizontal Deflection, Hollow Square Tube, [ALTERNATE DATA]')
xlim([0,2])
ylim([100,300])
hold on
df = linspace(100,300,100);
yeet = zeros(1,length(df)) + buck_s;
plot(delta,Ps,'LineWidth',1.3)
plot(yeet,df,'--g','LineWidth',1.3)
lb = zeros(1,length(0:.25:2)) + P_cr_square;
plot(0:.25:2,lb,'--k','LineWidth',1.3)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','SW')



figure
plot(dispRect2,voltRect2,'o','LineWidth',1.2)
title('Load vs. Horizontal Deflection, Rectangular Bar, [ALTERNATE DATA]')
hold on
ylim([120,160])
xlim([0,2])
df = linspace(120,160,100);
yeet = zeros(1,length(df)) + buck_r;
plot(delta,Pr,'LineWidth',1.3)
plot(yeet,df,'--g','LineWidth',1.3)
lb = zeros(1,length(0:.25:2)) + P_cr_rect;
plot(0:.25:2,lb,'--k','LineWidth',1.3)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force')
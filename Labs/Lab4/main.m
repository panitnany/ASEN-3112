%% ASEN 3112 Lab 4
% Written by : Panitnan Yuvanondha
% With modification from : (add your name if you add or modify anyhting)

% HouseKeeping
clc;clear;close all;

% Load .mat file of data
expData_bar2 = load('thin_bar_2.mat');
expData_bar2 = expData_bar2.data;
expData_bar1 = load('thin_bar_ 1.mat');
expData_bar1 = expData_bar1.data;
expData_rod2 = load('square_rod_2.mat');
expData_rod2 = expData_rod2.data;
expData_rod1 = load('square_rod_1.mat');
expData_rod1 = expData_rod1.data;

% First Column is Force but in mV
% Third Column is Deflection (in)


% Thin bar = Aluminum rectangular solid cross section
% Sqaure Rod = Aluminum square hollow cross section

%% Question 1 : Buckling Load
% Predict the buckling load for each

% Square Rod 
L_R = 11+3/8+0.4; % Length (in)
t = 0.0625; % Thickness (in)
l = 0.25; % outer dimension (length and width) (in)
E = 10*10^6; % Young's Modulus (psi)

I_R = (l^4)/12 + ((l-0.125)^4)/12; % moment of inerita of hollow rod

Pcr_R = ((pi^2*E*I_R)/L_R^2);

% Thin bar
L_b = 11+3/8+0.4;
l1 = 0.125;
l2 = 1;
E = 10*10^6;

I_b = 1/12*l1^3; % moment of inertia of the beam

Pcr_b = ((pi^2*E*I_b)/L_b^2);

% Experimentally determine the buckling load for each
 Bar_force = expData_bar2(:,2)/0.0023; % data conversion from V to lbs
 Rod_force = expData_rod2(:,2)/0.0023;
 
 Pcr_bar_exp = max(Bar_force);
 disp(Pcr_bar_exp)
 
 Pcr_rod_exp = max(Rod_force);
 disp(Pcr_rod_exp)
 
 %% Question 2 : Post-Buckling Behavior
 L=11+3/8+0.4;
 displace = linspace(0,2,10000); %Displacement from 0 to 2 in
v = @(del) del*sin((pi*(L/2)/L));

% Aluminum square hollow tube
E_sq = 10E6; %Young's Modulus [psi]
stress_sq = 35000; %Stress [psi]
t_sq = 0.0625; %Wall thickness [in]
I_sq = 0.0003051758; %Moment of inertia [in^4]
e_sq = stress_sq/E_sq; %Strain
Pcr_sq = (pi^2*E_sq*I_sq)/L^2; %Critical loads
k_sq = @(del) del*(pi^2/L^2)*sin((pi*(L/2)/L))*0.125;

% Aluminum rectangular solid
E_re = 10E6; %Young's Modulus [psi]
stress_re = 35000; %Stress [psi]
I_re = 1/12 * 0.125^3; %Moment of inertia [in^4]
e_re = stress_re/E_re; %Strain
Pcr_re = (pi^2*E_re*I_re)/L^2; %Critical loads
k_re = @(del) del*(pi^2/L^2)*sin((pi*(L/2)/L))*0.0625;

% Compute Buckling Loads for both Specimen
for i = 1:length(displace)
    v_x(i) = v(displace(i));
    k_sqx = k_sq(displace);
    k_rex = k_re(displace);
end
idx = 0.00001;

% Buckling loads
P_sq = displace(find(abs(k_sqx-double(e_sq))<idx,1)); 
P_re = displace(find(abs(k_rex-double(e_re))<idx,1));

% Analytical Prediction
displ = linspace(0,2,200); %displacement
P_square = Pcr_sq *( 1 + (pi^2 /(8*L^2))*displ.^2);
P_rec = Pcr_re*( 1 + (pi^2 /(8*L^2))*displ.^2);
 
figure
plot(expData_bar2(:,3),Bar_force,'o','LineWidth',1.3)
hold on
load_sq = linspace(100,300,100);
idx_sq = zeros(1,length(load_sq)) + P_sq;
plot(displ,P_square)
plot(idx_sq,load_sq,'--m','LineWidth',1.5)
xlim([0,2])
ylim([100,300])
title('Load vs Lateral Deflection - Tube')
xlabel('Deflection [in]')
ylabel('Load [lb]')

idx_sq2 = zeros(1,length(0:0.25:2)) + Pcr_sq;
plot(0:0.25:2,idx_sq2,'--k','LineWidth',1.5)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','E')

figure
plot(expData_rod2(:,3),Rod_force,'o','LineWidth',1.3)
title('Load vs Lateral Deflection - Beam')
xlim([0,2])
ylim([0,160])
hold on
load_re = linspace(0,160,100);
idx_re = zeros(1,length(load_re)) + P_re;
plot(displ,P_rec)
plot(idx_re,load_re,'--m','LineWidth',1.5)
idx_re2 = zeros(1,length(0:0.25:2)) + Pcr_re;
plot(0:0.25:2,idx_re2,'--k','LineWidth',1.5)
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','E')
xlabel('Deflection [in]')
ylabel('Load [lb]')
 
 
 %% Question 3 : Theoretical Design Study
 



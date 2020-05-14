%% House Keeping
clc
clear
close all

%% Constants
% Cross-section dimensions
L_tot = 4; %Total length of the entire truss from end to end [m]
d = 9.525*10^-3; %diameter of strutts (m)
t = 1.587*10^-3; %thickness of strutts (m)
P = 222.4; %total loads (N)
Lb = 250*10^-3; %Bay dimension, node location (m)
E = 68.9*10^9; % Modulus Elasticity of 6061-T6 Alluminum Alloy [Pa]
A = pi*(d^2)/4; % Area of the strutts [m^2]
Izz = 4*A*(Lb/2)^2; % Moment of Inertia of the strutts [m^4]
%% Read File
data = load('10lbincrements.txt');

% Sort data for each column and convert to SI units
L = data(:,1)*4.44822; %loading case [N]
F0 = data(:,2)*4.4482216; % [N]
F1 = data(:,3)*4.4482216; % [N]
F2 = data(:,4)*4.4482216; % [N]
F3D = data(:,5)*4.4482216; % [N]
LVDT = data(:,6)*25.4; % [mm]

% Average Data Points
datapt = 10; %data points per load
idx = 1;
for i = 1:(datapt)
    idx2 = i*10;
    L_ave(i) = mean(L(idx:idx2));
    F0_ave(i) = mean(F0(idx:idx2));
    F1_ave(i) = mean(F1(idx:idx2));
    F2_ave(i) = mean(F2(idx:idx2));
    F3D_ave(i) = mean(F3D(idx:idx2));
    LVDT_ave(i) = mean(LVDT(idx:idx2));
    idx = idx+datapt;
end

% Average of Loading and Unloading
for i = 1:(datapt)
    total_L(i) = (L_ave(i)+ L_ave((datapt)+1-i))/2;
    total_F0(i) = (F0_ave(i)+ F0_ave((datapt)+1-i))/2;
    total_F1(i) = (F1_ave(i)+ F1_ave((datapt)+1-i))/2;
    total_F2(i) = (F2_ave(i)+ F2_ave((datapt)+1-i))/2;
    total_F3D(i) = (F3D_ave(i)+ F3D_ave((datapt)+1-i))/2;
    total_LVDT(i) = (LVDT_ave(i)+ LVDT_ave((datapt)+1-i))/2;
end
%% Vertical Displacement
% Vertical displacement of the center of the truss (Equivalent Beam)
v = total_L.*L_tot^3/(48*E*Izz)*1000; %[mm]

%% Linearity of The Measurements vs External Load Magnitude
% Do with ANSYS
idx_a = 0.438095; %extrapolate ANSYS
a = idx_a/max(total_L);
val = total_L*a; %value for ANSYS
val_avg = mean(val);
%% Error
% Error for ANSYS
un_val = v-val_avg;
unc_val = sum(un_val.^2);
error_val = (unc_val/length(v))*100;

% Error for max displacement
un_disp = v-(val);
unc_disp = sum(un_disp.^2);
error_disp = (unc_disp/length(v))*100;

% Error for max displacement with ANSYS
un_disp2 = val-(v);
unc_disp2 = sum(un_disp2.^2);
error_disp2 = (unc_disp2/length(val))*100;
%% Plot

% Displacement
figure
hold on
plot(total_L,total_LVDT,'-o')
title('Vertical Displacement vs Applied Loading');
xlabel('Applied Load [N]');
ylabel('LVDT Displacement [mm]');

% Measured Force
figure
hold on
plot(total_L,total_F0,'-o')
plot(total_L,total_F1,'-o')
plot(total_L,total_F2,'-o')
plot(total_L,total_F3D,'-o')
xlabel('Applied Load [N]');
ylabel('Force [N]');
title('Force vs Applied Loads');
legend('F0','F1','F2','F3D','Location','southeast')

% Displacement
figure
hold on
plot(total_L,total_LVDT,'-o')
plot(total_L,v,'-o')
plot(total_L,val,'o-');
title('Vertical Displacement vs Applied Loading');
xlabel('Applied Load [N]');
ylabel('LVDT Displacement [mm]');
legend('Experimental Data','Equivalent Beam','ANSYS Measurements','Location','southeast');

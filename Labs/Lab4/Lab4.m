%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aufa Amirullah
% ASEN 3112 - Lab 4 (Problem 2)
% Instructor: Prof. Maute & Prof. Lopez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% House Keeping
clc
clear
close all

%% Read data
beam = load('lab4_data_group6-beam.mat');
tube = load('lab4_data_group6-square.mat');

beam = beam.data;
tube = tube.data;

% Data conversion from V to lbs
beam_Force = (beam(:,2)-beam(2,2))/0.0023;
tube_Force = (tube(:,2)-tube(2,2))/0.0023;

%% Constants
L = 10 + 5/16 + 0.4; %Length [in]
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

%% Compute Buckling Loads for both Specimen
for i = 1:length(displace)
    v_x(i) = v(displace(i));
    k_sqx = k_sq(displace);
    k_rex = k_re(displace);
end
idx = 0.00001;

% Buckling loads
P_sq = displace(find(abs(k_sqx-double(e_sq))<idx,1)); 
P_re = displace(find(abs(k_rex-double(e_re))<idx,1));

%% Analytical Prediction
displ = linspace(0,2,200); %displacement
P_square = Pcr_sq *( 1 + (pi^2 /(8*L^2))*displ.^2);
P_rec = Pcr_re*( 1 + (pi^2 /(8*L^2))*displ.^2);

%% Plot
figure
plot(tube(:,3),tube_Force,'o','LineWidth',1.3)
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
plot(beam(:,3),beam_Force,'o','LineWidth',1.3)
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
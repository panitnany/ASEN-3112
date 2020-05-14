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
expData_bar1(545:549,:)=[]; % delete bad data because deflection (column 3) shouldn't decrease 
% while the test's still operating
expData_rod2 = load('square_rod_2.mat');
expData_rod2 = expData_rod2.data;
expData_rod2 (610:613,:)= []; % same reason as bar 1
expData_rod1 = load('square_rod_1.mat');
expData_rod1 = expData_rod1.data;

% pull the first peice of data from each bin

%Bar 1
 counter = 1;
 for i = 0:1/16:2
    temp = find(expData_bar1(:,3) == i);
    mV_read(counter) = expData_bar1(temp(2),2)*1000;
    
    counter = counter+1;
     
 end
% 
% % convert to kg
lb_bar1 = mV_read/2.37;

%Bar 2
counter = 1;
 for i = 0:1/16:2
    temp = find(expData_bar2(:,3) == i);
    mV_read(counter) = expData_bar2(temp(2),2)*1000;
    
    counter = counter+1;
     
 end
% 
% % convert to kg
lb_bar2 = mV_read/2.37;

%Rod 1
counter = 1;
 for i = 0:1/16:2
    temp = find(expData_rod1(:,3) == i);
    mV_read(counter) = expData_rod1(temp(2),2)*1000;
    
    counter = counter+1;
     
 end
% 
% % convert to kg
lb_rod1 = mV_read/2.37;

% %Rod 2
% counter = 1;
%  for i = 0:1/16:2
%     temp = find(expData_rod2(:,3) == i);
%     mV_read(counter) = expData_rod2(temp(2),2)*1000;
%     
%     counter = counter+1;
%      
%  end
% 
% % convert to kg
lb_rod2 = mV_read/2.37;


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
 Bar_force2 = expData_bar2(:,2)/0.00237; % data conversion from V to lbs
 Rod_force2 = expData_rod2(:,2)/0.00237;
 
 Bar_force1 = expData_bar1(:,2)/0.00237; % data conversion from V to lbs
 Rod_force1 = expData_rod1(:,2)/0.00237;
 
 % Finding the maximum for all 4 data and then average them
 Pcr_bar_exp2 = max(Bar_force2);
 
 Pcr_rod_exp2 = max(Rod_force2);
 
 Pcr_bar_exp1 = max(Bar_force1);
 
 Pcr_rod_exp1 = max(Rod_force1);
 
 
 
 Pcr_rod_exp=mean([Pcr_rod_exp2 Pcr_rod_exp1]);
 
 Pcr_bar_exp=mean([Pcr_bar_exp2 Pcr_bar_exp1]);
 
 
 % Calculate the error
 RodError= abs((Pcr_rod_exp-Pcr_R)/Pcr_R)*100
 BarError= abs((Pcr_bar_exp-Pcr_b)/Pcr_b)*100
 
 %% Question 2 : Post-Buckling Behavior
 L=11+3/8+0.4;
displ = linspace(0,2,10000); %Displacement from 0 to 2 in
v = @(del) del*sin((pi*(L/2)/L)); 


 % Square Rod
 sigma_rod = 35000; % yield stress [psi]
 k_rod = @(del) del*(pi^2/L^2)*sin((pi*(L/2)/L))*0.125;
 
 % Thin bar
 sigma_bar= 35000; % yield stress [psi]
 e= sigma_bar/E; % Strain (Hooke's law)
 k_bar = @(del) del*(pi^2/L^2)*sin((pi*(L/2)/L))*0.0625;

 
 % Compute Buckling Loads for both Specimen
for i = 1:length(displ)
    v_x(i) = v(displ(i));
    k_barx = k_bar(displ);
    k_rodx = k_rod(displ);
end
idx = 0.00001;

% Buckling loads
P_b = displ(find(abs(k_barx-double(e))<idx,1)); 
P_r = displ(find(abs(k_rodx-double(e))<idx,1));

% Analytical Prediction
%displ = linspace(0,2,200); %displacement
P_rod = Pcr_R *( 1 + (pi^2 /(8*L^2))*displ.^2);
P_bar = Pcr_b*( 1 + (pi^2 /(8*L^2))*displ.^2);


% Plot (Rod)
figure
plot(expData_rod2(:,3),Rod_force2,'*'); hold on
plot(displ,P_rod,'LineWidth',4)
load_r = linspace(0,300,100);
dx_r = zeros(1,length(load_r)) + P_r;
plot(dx_r,load_r,'--m','LineWidth',1.5)
% dx_r is the answer to 2.1 (deflection prediction)
dx_r2 = zeros(1,length(0:0.25:2)) + Pcr_R;
plot(0:0.25:2,dx_r2,'--k','LineWidth',1.5)


xlabel('Horizontal Deflection (in)')
ylabel('Load (lbs)')
grid minor
title('Deflection VS Load (Square Rod)')
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','SouthEast')

% Plot (Bar)
figure
plot(expData_bar2(:,3),Bar_force2,'*'); hold on
plot(displ,P_bar,'LineWidth',4)
load_b = linspace(0,160,100);
dx_b = zeros(1,length(load_b)) + P_b;
plot(dx_b,load_b,'--m','LineWidth',1.5)
% dx_b is the answer to 2.1 (deflection prediction)
dx_b2 = zeros(1,length(0:0.25:2)) + Pcr_b;
plot(0:0.25:2,dx_b2,'--k','LineWidth',1.5)


xlabel('Horizontal Deflection (in)')
ylabel('Load (lbs)')
grid minor
title('Deflection VS Load (Thin Bar)')
legend('Experimental Data','Analytical Prediction','Predicted Start of Post-Buckling','Critical Force','Location','SouthEast')



 %% Question 3 : Theoretical Design Study
 
 



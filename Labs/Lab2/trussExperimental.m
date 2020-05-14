%% Truss Experimental Data Analysis
% 
% TO DO:
% - fix units (all metric?)
% - plots of reaction forces

% Housekeeping
clear; close all; clc

% load data for both weight increments
data5lb = load("ASEN 3112 - Sp20 - Lab 2 Data - 5 lb increments");
data10lb = load("ASEN 3112 - Sp20 - Lab 2 Data - 10 lb increments");

% sort data
load5 = data5lb(:,1); load10 = data10lb(:,1); %lbs
F0_5 = data5lb(:,2); F0_10 = data10lb(:,2); %lbf
F1_5 = data5lb(:,3); F1_10 = data10lb(:,3); %lbf
F2_5 = data5lb(:,4); F2_10 = data10lb(:,4); %lbf
F3D_5 = data5lb(:,5); F3D_10 = data10lb(:,5); %lbf
LVDT5 = data5lb(:,6); LVDT10 = data10lb(:,6); %in

% conversions
% 50 lbs = 222.4 N, use to convert lbs -> N 
load5 = load5*222.4/50; load10 = load10*222.4/50;

% linear regression
[p5,S5] = polyfit(load5,LVDT5,1);
xFit = [0:5:50]*222.4/50; %cut to ascending half of load data
[yFit5,error5] = polyval(p5, xFit, S5);
maxError5 = max(error5) %max and average error values
meanError5 = mean(error5)

[p10,S10] = polyfit(load10,LVDT10,1);
xFit = [0:5:50]*222.4/50; %cut to ascending half of load data
[yFit10,error10] = polyval(p10, xFit, S10);
maxError10 = max(error10) %max and average error values
meanError10 = mean(error10)

% plot displacement over time
figure
hold on
plot(load5,LVDT5,'o', xFit,yFit5)
xlabel('Load [N]'); ylabel('Linear displacement [in]');
title('Displacement vs. Load for 5 lb. Increments');

figure
hold on
plot(load10,LVDT10,'o', xFit,yFit10)
xlabel('Load [N]'); ylabel('Linear displacement [in]');
title('Displacement vs. Load for 10 lb. Increments');

%% ASEN 3112 Lab 1
% Section 2.2
% Matthew Bridges (105070470)

% Housekeeping
clear;
clc;
close all;

%% Loading Data and Organizing data

ctwdata= importfile('CTW.txt');
otwdata= importfile('OTW.txt');

time_otw = otwdata(:,1);
time_ctw = ctwdata(:,1);

% Organzie the data
% Closed thin walled
gammactw= deg2rad(ctwdata(:,2)); % shear strain
phictw= ctwdata(:,3)*pi/180; % total twist angle in rad
Tctw= ctwdata(:,4); % applied torque in-lb

% Open thin walled
% machine twist angle reset for zero at first data point (rad)
strain_otw= (otwdata(:,2)-otwdata(1,2))*pi/180; 
% extensometer twist angle (rad)
strain_ext_otw= (otwdata(:,3)-otwdata(1,3))*pi/180; 
Totw= otwdata(:,4); % applied torque in-lb
dx = otwdata(:,5);

%% Constant Variables
G=3.75*10^6; % psi
Re= 3/8; % in. for exterior radius
t= 1/16; % in. for thickness
Ri= Re-t; % in. for inner radius
Rav = (Ri+Re)/2; % in. for average radius 
L= 9; % in. length of the tube
Le= 1; % in. length of extensometer
alpha= 1/3; 

%% Section 2.2: OTW

% Section 2.2.1
% Calculating twist angle using shear strain
phi_otw = strain_otw*(L/Re);
phi_ext_otw = strain_ext_otw*(L/Re);

% Plot
figure(1)
hold on
xlabel('Torque (lb-in)')
ylabel('Shear Strain')
title('Torque VS Shear strain (OTW)')
plot(Totw,strain_otw)
grid minor
plot(Totw,strain_ext_otw)
legend('Machine calculated shear strain','Extensometer calculated shear strain')

%% Finding torsional rigidty
% Preparing variables for least fit approx
strain_otw = strain_otw(1000:end);
Totw = Totw(1000:end);
strain_ext_otw = strain_ext_otw(1000:end);
% Least squares approxs
p_otw = polyfit(strain_otw,Totw*Re,1);
p_ext_otw = polyfit(strain_ext_otw,Totw*Re,1);

% Plotting Strain and Torque with best fit line
figure(2)
hold on
scatter(strain_otw,Totw*Re,'*b')
plot(strain_otw,(p_otw(1)*strain_otw+p_otw(2)),'r')
grid minor
xlabel('Strain')
ylabel('Torque')
title('Torque VS Shear strain (OTW)')
legend('Machine data','Best fit')

GJ_otw = -p_otw(1) % Torsional rigidity is the slope of T/phi
std_otw = std(Totw./strain_otw); % standard deviation of GJ
sigma_GJotw = std_otw/sqrt(length(Totw)) % error in GJ

% Repeatting for extensometer data
figure(3)
hold on
scatter(strain_ext_otw,Totw*Re,'*b')
plot(strain_ext_otw,(p_ext_otw(1)*strain_ext_otw+p_otw(2)),'r')
grid minor
xlabel('Strain')
ylabel('Torque')
title('Torque VS Shear strain (OTW)')
legend('Extensometer data','Best fit')

GJ_ext_otw = -p_ext_otw(1)
std_ext = std(Totw./strain_ext_otw);
sigma_GJotw_ext = std_ext/sqrt(length(Totw))



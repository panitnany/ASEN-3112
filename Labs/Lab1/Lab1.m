%% House Keeping
clc
clear; close all

%% Constants
Re = 3/8; %outer radius (in)
t = 1/16; %wall thickness (in)
Ri = Re - t; %in
t = 1/16; %wall thickness (in)
G = 3.75 *10^6; %shear modulus
Le = 1; %extensometer gauge length (in)
L = 9;
alpha = 1/3;

%% Load Files
CTW = importfile('CTW.txt');
OTW = importfile('OTW.txt');
%% Analysis of the Close Thin Wall Specimen
%% Parse Data
phiCTW = (CTW(:,2) - CTW(1,2))*pi/180;
TorqueCTW = CTW(:,4);

[~,maxStrainExt] = max(abs(CTW(:,3)));
[~,maxStrainPhi] = max(abs(phiCTW));

StrainCTWExt = CTW(1:maxStrainExt,3)*pi/180;
StrainCTWphi = phiCTW(1:maxStrainPhi)*(Re-(t/2))/L;

%% Least Squares Estimation
xCTW = polyfit(TorqueCTW(1:maxStrainExt),StrainCTWExt,1);
xCTW2 = polyfit(TorqueCTW(1:maxStrainPhi),StrainCTWphi,1);

%% Plotting

plot(TorqueCTW(1:length(StrainCTWExt)),-StrainCTWExt,'b')
hold on 
plot(TorqueCTW(1:length(StrainCTWphi)),-StrainCTWphi,'r')
title('Torque vs. Strain from Extensometer and Twist Angle')
xlabel('Torque (Lb-in)')
ylabel('Strain (in/in)')
plot(TorqueCTW, -xCTW(1)*TorqueCTW + xCTW(2),'g--')
plot(TorqueCTW, -xCTW2(1)*TorqueCTW + xCTW2(2),'m--')
legend('Strain from Extensometer', 'Strain from Twist Angle','Extensometer Best Fit Line',...
    'Twist Angle Best Fit Line')

%% Torsional Rigidity
TRext = Re/-xCTW(1); %Lb*in^3
TRphi = Re/-xCTW2(1); %Lb*in^3

%% Analysis of the Open Thin Wall Specimen

% 2.2.1 - Plot Torque vs Strain From The Extensometer and Testing Machine
% Define data from the extensometer
t_o = OTW(:,1); %time [s]
phi_ext= deg2rad(OTW(:,2)); %twist angle [rad]
phi_ext = phi_ext-phi_ext(1);
gamma_ext = deg2rad(OTW(:,3)); %shear strain [rad]
T_ext = OTW(:,4); %Torque [lb-in]

% Best estimate i.e. Theoretical Prediction
gamma_the = (phi_ext.*t)/Le;
x_the = polyfit(T_ext,gamma_the,1);
GJ_the = t/x_the(1); %Torsional rigidity

%% 2.2.2 - Least Fitting Square to calculate the Torsional Rigidity
% From Extensometer
x_ext = polyfit(T_ext,gamma_ext,1);
GJ_ext = t/x_ext(1);
% Calculating the uncertainty
std_devext = std(gamma_ext);
unc_ext = std_devext/sqrt(numel(gamma_ext));

% From Testing Machine
gamma_tes = (phi_ext.*t)/L; %shear strain
x_tes = polyfit(T_ext,gamma_tes,1);
GJ_tes = t/x_tes(1);

% Calculating the uncertainty
std_devtes = std(gamma_tes);
unc_tes = std_devtes/sqrt(numel(gamma_tes));

%% Plot Torques vs Shear Strain
figure
hold on

%Theoretical Models
plot(T_ext,gamma_ext); %by the extensometer
plot(T_ext,gamma_tes); %by the total rotation of testing machine

%Experimental Model
plot(T_ext,(x_ext(1)*T_ext+x_ext(2))); %by the extensometer
plot(T_ext,(x_tes(1)*T_ext+x_tes(2))); %by the total rotation of testing machine

xlabel('Torque [lb-in]');
ylabel('Shear strain [rad]');
title('Torque vs Shear Strain (OTW)');
legend('Experimental Model 1 - by The Extensometer','Experimental Model 2 - by The Total Rotation','Theoretical Model 1 - by The Extensometer','Theoretical Model 2 - by The Total Rotation')

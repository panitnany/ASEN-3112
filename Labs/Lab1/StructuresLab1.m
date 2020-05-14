% housekeeping 
clc;clear;close all

%% Loading Data and Organizing data

ctwdata= importfile('CTW.txt');
otwdata= importfile('OTW.txt');

% Organzie the data
% Closed thin walled
gammactw= (ctwdata(:,2))*pi/180; % shear strain
phictw= (ctwdata(:,3)*pi/180); % total twist angle in rad
Tctw= ctwdata(:,4); % applied torque in-lb

% Open thin walled     
gammaotw= otwdata(:,2)*pi/180; % shear strain
phiotw= otwdata(:,3)*pi/180; % total twist angle in rad
Totw= otwdata(:,4); % applied torque in-lb


%% Constant Variables
G=3.75*10^6; % psi
Re= 3/8; % in. for exterior radius
t= 1/16; % in. for thickness
Ri= Re-t; %in. for inner radius
L= 9; % in. length of the tube
Le= 1; % in. length of extensometer
alpha= 1/3; 

%% Section 2.1 : CTW 

% Section 2.1.1
% Calculating strain by using total rotation
strain = (phictw(:,1)*Re)/L;

% Plot
figure(1)
hold on
xlabel('Torque (lb-in)')
ylabel('Shear Strain')
title('Torque VS Shear strain (CTW)')
plot(Tctw,strain,'r') % by total rotation (Calculation)
grid minor
plot(Tctw,phictw,'k') % by extensomter (provided)


% Section 2.1.2 polyfit (?????

lsgamma= polyfit(Tctw,phictw,1); % least square fitting for shear strain provided
lsstrain = polyfit(Tctw,strain,1); % least sqaure fitting for shear strain calculated

hold on
plot(Tctw, lsgamma(1)*Tctw+lsgamma(2),'m--', 'LineWidth',2)
plot(Tctw, lsstrain(1)*Tctw+ lsstrain(2),'b--','LineWidth',2)
legend('Calculated shear strain (from twist angle)','shear strain provided'...
    ,'Best fit line for provided value','Best fit line for twist angle','Location','southwest')


%Section 2.1.3 
% Calculate J CTW (theoretical)
R=(Re+Ri)/2; % midline
Ae= pi*R^2; % Enclosed Area
Jctw= 4*Ae^2/(2*pi*R/t);

% Calculate J exact theory (theoretical)
Jex= pi/2*(Re^4-Ri^4);

% Calculate torsion regidity (theoretical)
TRctw= G*Jctw; % closed thoery
TRex= G* Jex; % exact theory

% Calculate torsion regidity (experimental)
TRctwexp= -R/lsgamma(1); % from shear strain provided

% Calculate torsion regidity (experimental)
TRctwexp1= -R/lsstrain(1); % from shear strain calculated


%% Section 2.2: OTW

% Section 2.2.1
% Calculating strain by using total rotation
strainotw = (phiotw(:,1)*t)/L;

% Plot
figure(2)
hold on
xlabel('Torque (lb-in)')
ylabel('Shear Strain')
title('Torque VS Shear strain (OTW)')
plot(Totw,strainotw)
grid minor
plot(Totw,phiotw)
ylim([-2.5*10^-3,1*10^-3])


% Section 2.2.2

lsgammaotw= polyfit(Totw,phiotw,1); % least square fitting for shear strain provided
lsstrainotw = polyfit(Totw,strainotw,1); % least sqaure fitting for shear strain calculated

hold on
plot(Totw, lsgammaotw(1)*Totw+lsgammaotw(2),'g--')
plot(Totw, lsstrainotw(1)*Totw+ lsstrainotw(2),'c--')
legend('Calculated shear strain (from twist angle)','shear strain provided'...
    ,'Best fit line for provided value','Best fit line for twist angle')


%Section 2.1.3 
% Calculate J OTW (theoretical)
Jotw= alpha*2*pi*R*t^3;

% Calculate torsion regidity (theoretical)
TRotw= G*Jotw;

% Calculate torsion regidity (experimental)
TRotwexp= R/lsgammaotw(1); % from shear strain provided
TRotwexp1= R/lsstrainotw(1); % from shear strain calculated


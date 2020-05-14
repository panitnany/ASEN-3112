%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nate Lee
% ASEN 3112
% Lab 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Good Houskeeping
clear all 
close all

%% Calculations

% Set L Vector
Lvec = [10:0.05:12];

Eflat = 10000000; %psi
Esquare = 10000000; %psi

Iflat = 0.0001627604; %in^4
Isquare = 0.0003051758; %in^4

% Simply Supported
PsquareSS = (pi^2*Esquare*Isquare)./(Lvec).^2;
PflatSS = (pi^2*Eflat*Iflat)./(Lvec).^2;

% Fixed Fixed
PsquareFF = (pi^2*Esquare*Isquare)./(Lvec.*0.5).^2;
PflatFF = (pi^2*Eflat*Iflat)./(Lvec.*0.5).^2;

%% Plotting Section

% Fixed
figure
hold on
plot(Lvec,PsquareFF,'linewidth',2)
plot(Lvec,PflatFF,'linewidth',2)
grid on
legend('Square Beam','Rectangular Beam')
xlabel('Length (in)')
ylabel('Load (lb)')
title('Critical Load with Respect to Length for a Fixed-Fixed Beam')

%figure
figure
hold on
plot(Lvec,PsquareSS,'linewidth',2)
plot(Lvec,PflatSS,'linewidth',2)
grid on
legend('Square Beam','Rectangular Beam')
xlabel('Length (in)')
ylabel('Load (lb)')
title('Critical Load with Respect to Length for a SS-SS Beam')


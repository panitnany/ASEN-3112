clc; clear; close all;
%% Q1
% predicting Pcr and comparing to data
beam = load('lab4_data_group6-beam.mat');
beam = beam.data;
tube = load('lab4_data_group6-square.mat');
tube = tube.data;

% 2.3 mV per pound

Beam_force = (beam(:,2)-beam(2,2))/0.0023;
Tube_force = (tube(:,2)-tube(2,2))/0.0023;


%% square tube theoretical

L = 10+5/16+0.4;
t = 0.0625;
l = 0.25;
E = 10*10^6;

I_t = (l^4)/12 + ((l-0.125)^4)/12; % moment of inerita of the tube

Pcr_t = ((pi^2*E*I_t)/L^2);

%% beam theoretical

L_b = 11 + 7/6 + 0.4;
l1 = 0.125;
l2 = 1;
E = 10*10^6;

I_b = 1/12 * 0.125^3; % moment of inertia of the beam

Pcr_b = ((pi^2*E*I_b)/L_b^2);

%% Finding Pcr from the experimentla data

Pcr_tube_exp = max(Tube_force);
disp(Pcr_tube_exp)

% the Pcr value is not clear therefor an average is taken over the region
Pcr_beam_exp = mean(Beam_force(24:85));
disp(Pcr_beam_exp)







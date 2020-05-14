% this code analyzes the LA data for the 4th Lab of ASEN 3112
% William Butler
% 3/30/20

%% Housekeeping
clear all
close all
clc

%% Import data
load('square_rod_1.mat')
sq1_data = data;
load('square_rod_2.mat')
sq2_data = data;
load('thin_bar_ 1.mat')
bar1_data = data;
load('thin_bar_2.mat')
bar2_data = data;

%% Knowns
% Rectangular rod material properties
L_bar = 11 + 3/8;
% Square rod material properties
L_sq = 11 + 3/8;

%% Analyze Sqaure 1 Data
figure
plot(sq1_data(:,3),sq1_data(:,2),'.')
xlabel('Displacement Bin [in.]')
ylabel('Voltage [mV]')
title('Hollow Square Rod 1: Raw Data (Voltage vs. Displacement)')

% pull the first peice of data from each bin
counter = 1;
for i = 0:1/16:2
   temp = find(sq1_data(:,3) == i);
   mV_read(counter) = sq1_data(temp(2),2)*1000;
   
   counter = counter+1;
    
end

% convert to kg
lb_sq1 = mV_read/2.37;
%convert to newtons
% N_sq1 = kg_sq1*9.8;

% generate a fit structure
sq1_fit = fit([0:1/16:2]',lb_sq1', 'smoothingspline');

% find the slope of the elastic region
%sq1_el_fit = polyfit([0:1/16:2/16]',lb_sq1(1:3)', 1);

% generate yeild inclined line
%sq1_el_disp = linspace(0,3/16,1000);
%sq1_el_line = (sq1_el_fit(1))*(sq1_el_disp-(L_sq*.002));

% find critical point
%[sq1_F_buck, I2] = max(lb_sq1);
%sq1_disp_buck = (1/16)*I2;

figure
hold on
% plot([0:1/16:2],N_sq1,'o')
plot(sq1_fit, [0:1/16:2],lb_sq1, 'bo')
%plot(sq1_el_disp, sq1_el_line, 'm--')
%line(sq1_disp_buck, 'k--');
title('Force vs. Displacement for Hollow Square Rod 1')
ylabel('Force [lb]')
xlabel('Displacement [in.]')
%legend('Raw Data', 'Fitted Line', 'Elastic Yield Line', 'Critical Load')
legend('Raw Data', 'Fitted Line')
ylim([0,230])
grid on

% find the intersection of the inclined yield line
%I = find(sq1_el_line>=199.3);
%sq1_F_yield = sq1_el_line(I(1));
%sq1_disp_yield = sq1_el_disp(I(1));



% generate outputs
%fprintf('******************************************************************\n')
%fprintf('Square Rod 1 Experimental Outputs:\n')
%fprintf('The Yield Force is: ~%0.3f lbF\n',sq1_F_yield)
%fprintf('The Yield Displacement is: ~%0.3f in\n', sq1_disp_yield)
%fprintf('The Critical Load is: ~%0.3f lbF\n',sq1_F_buck)
%fprintf('The Critical Displacement is: ~%0.3f in\n',sq1_disp_buck)
%fprintf('******************************************************************\n')
%% Analyze Sqaure 2 Data
figure
plot(sq2_data(:,3),sq2_data(:,2),'.')
xlabel('Displacement Bin [in.]')
ylabel('Voltage [mV]')
title('Hollow Square Rod 2: Raw Data (Voltage vs. Displacement)')


% pull the first peice of data from each bin
counter = 1;
for i = 0:1/16:2
   try
   temp = find(sq2_data(:,3) == i);
   mV_read(counter) = sq2_data(temp(2),2)*1000;
   
   counter = counter+1;
   
   catch
       fprintf('Square 2 Data:\n')
       fprintf('Displacement Bin Missing, Skipping Bin: %0.4f\n',i)
   end
end

% convert to kg
lb_sq2 = mV_read/2.37;
%convert to newtons
% N_sq2 = kg_sq2*9.8;

% generate a fit structure
sq2_fit = fit([0:1/16:2]',lb_sq2', 'smoothingspline');

% find the slope of the elastic region
%sq2_el_fit = polyfit([0:1/16:2/16]',lb_sq2(1:3)', 1);

% generate yeild inclined line
%sq2_el_disp = linspace(0,3/16,1000);
%sq2_el_line = (sq2_el_fit(1))*(sq2_el_disp-(L_sq*.002));

% find critical point
%[sq2_F_buck, I2] = max(lb_sq2);
%sq2_disp_buck = (1/16)*I2;

figure
hold on
% plot([0:1/16:2],N_sq2,'o')
plot(sq2_fit, [0:1/16:2],lb_sq2, 'bo')
%plot(sq1_el_disp, sq1_el_line, 'm--')
%line(sq1_disp_buck, 'k--');
title('Force vs. Displacement for Hollow Square Rod 2')
ylabel('Force [lb]')
xlabel('Displacement [in.]')
%legend('Raw Data', 'Fitted Line', 'Elastic Yield Line', 'Critical Load')
legend('Raw Data', 'Fitted Line')
ylim([0,230])
grid on

% find the intersection of the inclined yield line
%I = find(sq2_el_line>=211.8);
%sq2_F_yield = sq2_el_line(I(1));
%sq2_disp_yield = sq2_el_disp(I(1));

% generate outputs
%fprintf('******************************************************************\n')
%fprintf('Square Rod 2 Experimental Outputs:\n')
%fprintf('The Yield Force is: ~%0.3f lbF\n',sq2_F_yield)
%fprintf('The Yield Displacement is: ~%0.3f in\n', sq2_disp_yield)
%fprintf('The Critical Load is: ~%0.3f lbF\n',sq2_F_buck)
%fprintf('The Critical Displacement is: ~%0.3f in\n',sq2_disp_buck)
%fprintf('******************************************************************\n')

%% Analyze Bar 1 Data
figure
plot(bar1_data(:,3),bar1_data(:,2),'.')
xlabel('Displacement Bin [in.]')
ylabel('Voltage [mV]')
title('Rectangular Bar 1: Raw Data (Voltage vs. Displacement)')

% pull the first peice of data from each bin
counter = 1;
for i = 0:1/16:2
   try
   temp = find(bar1_data(:,3) == i);
   mV_read(counter) = bar1_data(temp(2),2)*1000;
   
   counter = counter+1;
   
   catch
       fprintf('Bar 1 Data:\n')
       fprintf('Displacement Bin Missing, Skipping Bin: %0.4f\n',i)
   end
end

% convert to kg
lb_bar1 = mV_read/2.37;
%convert to newtons
% N_bar1 = kg_bar1*9.8;

% generate a fit structure
bar1_fit = fit([0:1/16:2]',lb_bar1', 'smoothingspline');

% find the slope of the elastic region
%bar1_el_fit = polyfit([0:1/16:2/16]',lb_bar1(1:3)', 1);

% generate yeild inclined line
%bar1_el_disp = linspace(0,6/16,1000);
%bar1_el_line = (bar1_el_fit(1))*(bar1_el_disp-(L_bar*.002));

% find critical point
%[bar1_F_buck, I2] = max(lb_bar1);
%bar1_disp_buck = (1/16)*I2;

figure
hold on
% plot([0:1/16:2],N_bar1,'o')
plot(bar1_fit, [0:1/16:2],lb_bar1, 'bo')
%plot(bar1_el_disp, bar1_el_line, 'm--')
%line(bar1_disp_buck, 'k--');
title('Force vs. Displacement for Rectangular Bar 1')
ylabel('Force [lb]')
xlabel('Displacement [in.]')
%legend('Raw Data', 'Fitted Line', 'Elastic Yield Line', 'Critical Load')
legend('Raw Data', 'Fitted Line')
ylim([0,120])
grid on

% find the intersection of the inclined yield line
%I = find(bar1_el_line>=90);
%bar1_F_yield = bar1_el_line(I(1));
%bar1_disp_yield = bar1_el_disp(I(1));

% generate outputs
%fprintf('******************************************************************\n')
%fprintf('Rectangular Rod 1 Experimental Outputs:\n')
%fprintf('The Yield Force is: ~%0.3f lbF\n',bar1_F_yield)
%fprintf('The Yield Displacement is: ~%0.3f in\n', bar1_disp_yield)
%fprintf('The Critical Load is: ~%0.3f lbF\n',bar1_F_buck)
%fprintf('The Critical Displacement is: ~%0.3f in\n',bar1_disp_buck)
%fprintf('******************************************************************\n')

%% Analyze Bar 2 Data
figure
plot(bar2_data(:,3),bar2_data(:,2),'.')
xlabel('Displacement Bin [in.]')
ylabel('Voltage [mV]')
title('Rectangular Bar 2: Raw Data (Voltage vs. Displacement)')

% pull the first peice of data from each bin
counter = 1;
for i = 0:1/16:2
   try
   temp = find(bar2_data(:,3) == i);
   mV_read(counter) = bar2_data(temp(2),2)*1000;
   
   counter = counter+1;
   
   catch
       fprintf('Bar 2 Data:\n')
       fprintf('Displacement Bin Missing, Skipping Bin: %0.4f\n',i)
   end
end

% convert to kg
lb_bar2 = mV_read/2.37;
%convert to newtons
% N_bar2 = kg_bar2*9.8;

% generate a fit structure
bar2_fit = fit([0:1/16:2]',lb_bar2', 'smoothingspline');

% find the slope of the elastic region
%bar2_el_fit = polyfit([0:1/16:2/16]',lb_bar2(1:3)', 1);

% generate yeild inclined line
%bar2_el_disp = linspace(0,6/16,1000);
%bar2_el_line = (bar2_el_fit(1))*(bar2_el_disp-(L_bar*.002));

% find critical point
%[bar2_F_buck, I2] = max(lb_bar2);
%bar2_disp_buck = (1/16)*I2;

figure
hold on
% plot([0:1/16:2],N_bar2,'o')
plot(bar2_fit, [0:1/16:2],lb_bar2, 'bo')
%plot(bar1_el_disp, bar1_el_line, 'm--')
%line(bar1_disp_buck, 'k--');
title('Force vs. Displacement for Rectangular Bar 2')
ylabel('Force [lb]')
xlabel('Displacement [in.]')
%legend('Raw Data', 'Fitted Line', 'Elastic Yield Line', 'Critical Load')
legend('Raw Data', 'Fitted Line')
ylim([0,120])
grid on

% find the intersection of the inclined yield line
%I = find(bar2_el_line>=95.23);
%bar2_F_yield = bar2_el_line(I(1));
%bar2_disp_yield = bar2_el_disp(I(1));

% generate outputs
%fprintf('******************************************************************\n')
%fprintf('Rectangular Rod 2 Experimental Outputs:\n')
%fprintf('The Yield Force is: ~%0.3f lbF\n',bar2_F_yield)
%fprintf('The Yield Displacement is: ~%0.3f in\n', bar2_disp_yield)
%fprintf('The Critical Load is: ~%0.3f lbF\n',bar2_F_buck)
%fprintf('The Critical Displacement is: ~%0.3f in\n',bar2_disp_buck)
%fprintf('******************************************************************\n')

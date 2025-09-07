clear all; close all; clc;

ricker_data = dlmread('1_Ricker_model.txt', '', 1, 0);
generation = ricker_data(:,1);
x_values = ricker_data(:,2:5);

% Plot Ricker model population
figure(1);
plot(generation, x_values(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'r=1.9');
hold on;
plot(generation, x_values(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'r=2.9');
plot(generation, x_values(:,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'r=3.3');
plot(generation, x_values(:,4), 'm-', 'LineWidth', 1.5, 'DisplayName', 'r=3.6');
hold off;

xlabel('Generation');
ylabel('Fractional Population (x)');
title('Ricker Model Population Dynamics');
legend('Location', 'best');


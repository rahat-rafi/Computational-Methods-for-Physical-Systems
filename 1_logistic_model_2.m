clear all; close all; clc;

logistic_data = dlmread('1_logistic_model_2.txt', '', 1, 0);
generation = logistic_data(:,1);
x_values = logistic_data(:,2:5);
p_values = logistic_data(:,6:9);
K_values = logistic_data(:,10);

% Plot fractional population (x)
figure(1);
plot(generation, x_values(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'r=1.9');
hold on;
plot(generation, x_values(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'r=2.9');
plot(generation, x_values(:,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'r=3.3');
plot(generation, x_values(:,4), 'm-', 'LineWidth', 1.5, 'DisplayName', 'r=3.6');
plot(generation, K_values, 'k--', 'LineWidth', 2, 'DisplayName', 'Carrying Capacity K');
hold off;

xlabel('Generation');
ylabel('Fractional Population (x)');
title('Logistic Model with Decreasing K - Fractional Population');
legend('Location', 'best');


figure(2);
plot(generation, p_values(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'r=1.9');
hold on;
plot(generation, p_values(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'r=2.9');
plot(generation, p_values(:,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'r=3.3');
plot(generation, p_values(:,4), 'm-', 'LineWidth', 1.5, 'DisplayName', 'r=3.6');
plot(generation, K_values, 'k--', 'LineWidth', 2, 'DisplayName', 'Carrying Capacity K');
hold off;

xlabel('Generation');
ylabel('Absolute Population (p)');
title('Logistic Model with Decreasing K - Absolute Population');
legend('Location', 'best');

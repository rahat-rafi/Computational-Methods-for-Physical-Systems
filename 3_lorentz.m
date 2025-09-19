clear; clc;

fname = '3_lorentz_dt_1.00000.dat';

data = dlmread(fname, '', 1, 0);

time  = data(:,1);
xpos  = data(:,2);
ypos  = data(:,3);
zpos  = data(:,4);
err_x = data(:,6);

base = regexprep(fname, '\.dat$', '');

% -------- Plot 1: y vs x --------
figure(1); clf;
plot(xpos, ypos, 'b-');
axis equal; grid on;
xlabel('x', 'FontWeight', 'bold');
ylabel('y', 'FontWeight', 'bold');
title('Trajectory: y vs x', 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5);
print([base '_y_vs_x.png'], '-dpng', '-r300');

% -------- Plot 2: z vs time --------
figure(2); clf;
plot(time, zpos, 'r-');
grid on;
xlabel('time', 'FontWeight', 'bold');
ylabel('z', 'FontWeight', 'bold');
title('z vs time', 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5);
print([base '_z_vs_time.png'], '-dpng', '-r300');

% -------- Plot 3: log10(|Δx|) vs time --------
figure(3); clf;
plot(time, log10(abs(err_x)), 'k-');
grid on;
xlabel('time', 'FontWeight', 'bold');
ylabel('log_{10}(|\Delta x|)', 'FontWeight', 'bold');
title('Error growth in x', 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5);
print([base '_log10_errx_vs_time.png'], '-dpng', '-r300');


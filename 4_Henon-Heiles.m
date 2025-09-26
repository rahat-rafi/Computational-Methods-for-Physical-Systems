clear; clc;

% Load files (skip header line "# ...")
h1 = dlmread("4_henon_0.0900.dat", "", 1, 0);
h2 = dlmread("4_henon_0.1250.dat", "", 1, 0);
h3 = dlmread("4_henon_0.1482.dat", "", 1, 0);

% ---------- Trajectory Plot (y vs x) ----------
fig1 = figure(1); clf;
plot(h3(:,2), h3(:,3), 'r-', 'LineWidth',1.0, 'DisplayName','T0=0.1482'); hold on;
plot(h2(:,2), h2(:,3), 'g-', 'LineWidth',1.0, 'DisplayName','T0=0.1250');
plot(h1(:,2), h1(:,3), 'b-', 'LineWidth',1.0, 'DisplayName','T0=0.0900');

xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("Hénon–Heiles trajectories (y vs x)", "FontSize", 16, "FontWeight", "bold");
legend("location","northwest", "FontSize", 12);
axis equal; grid on;
set(gca, "FontSize", 12, "LineWidth", 1);
print(fig1, '4_henon_y_vs_x.png', '-dpng', '-r300');

% ---------- Phase-space Plot (dy/dt vs y) ----------
fig2 = figure(2); clf;
plot(h3(:,3), h3(:,5), 'r.', 'MarkerSize',8, 'DisplayName','T0=0.1482'); hold on;
plot(h2(:,3), h2(:,5), 'g.', 'MarkerSize',8, 'DisplayName','T0=0.1250');
plot(h1(:,3), h1(:,5), 'b.', 'MarkerSize',8, 'DisplayName','T0=0.0900');

xlabel("y", "FontSize", 14);
ylabel("dy/dt (p_y)", "FontSize", 14);
title("Phase-space plot: dy/dt vs y", "FontSize", 16, "FontWeight", "bold");
legend("location","northwest", "FontSize", 12);
grid on;
set(gca, "FontSize", 12, "LineWidth", 1);
print(fig2, '4_henon_Py_vs_y.png', '-dpng', '-r300');


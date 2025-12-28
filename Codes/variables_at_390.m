clear all; close all; clc;

% Load snapshot file:
% 1: iter   2: j-index   3: y   4: u   5: v   6: psi   7: vort
x1 = load('case2_psi_ix390_snap_w100.dat');

% Iterations to plot
iters = [20000 40000 60000 80000 100000 ...
         120000 140000 160000 180000 200000];

% Explicit RGB colors (one row per iteration)
% 20000  -> blue
% 40000  -> red
% 60000  -> dark green
% 80000  -> magenta
% 100000 -> cyan
% 120000 -> black
% 140000 -> yellow
% 160000 -> gray
% 180000 -> brownish orange
% 200000 -> dark purple
colors = [
    0.0   0.0   1.0 ;   % blue
    1.0   0.0   0.0 ;   % red
    0.0   0.5   0.0 ;   % dark green
    1.0   0.0   1.0 ;   % magenta
    0.0   1.0   1.0 ;   % cyan
    0.0   0.0   0.0 ;   % black
    1.0   1.0   0.0 ;   % yellow
    0.5   0.5   0.5 ;   % gray
    0.8   0.4   0.0 ;   % brownish orange
    0.3   0.0   0.5     % dark purple
];

figure; hold on;

for k = 1:length(iters)
    it = iters(k);

    % Rows for this iteration
    idx = (x1(:,1) == it);

    y   = x1(idx, 2);   % physical y
    psi = x1(idx, 6);   % stream function ψ at ix = 390

    plot(y, psi, 'LineWidth', 1.5, 'Color', colors(k,:));
end

hold off;
% LaTeX-style axis labels
xlabel('y',  'FontSize', 16);
ylabel('\psi', 'FontSize', 16);

xlim([0 450]);
ylim([0 450]);

set(gca, 'FontSize', 16, 'Box', 'on');

saveas(gcf, 'Bc2_psi_vs_y_w100', 'png');


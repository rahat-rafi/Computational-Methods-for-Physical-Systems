clear all; close all;

data = load('extra1_wswitch_field_200K_step10.dat');

x = data(:,1);
y = data(:,2);
u = data(:,3);
v = data(:,4);

figure(1)

% velocity field (red arrows)
quiver(x, y, u, v, 0.5, 'r');   % 0.5 is the scale factor – adjust if needed
hold on

% --- rectangular obstruction: white interior, blue dashed border ---
% In your grid: block is i = 401..451, j = 1..51
% and you wrote x = i-1, y = j-1 in the Fortran code,
% so the physical extents are x = 400..450, y = 0..50.

xb = 400;      % left x
yb = 1;        % bottom y
wb = 50;       % width  (x direction)
hb = 50;       % height (y direction)

rectangle('Position', [xb yb wb hb], ...
          'EdgeColor', 'b', ...
          'LineStyle', ':', ...   % dotted / dashed blue outline
          'LineWidth', 1.5, ...
          'FaceColor', 'w');      % white fill to hide arrows inside

axis([0 855 1 451])
%set(gca, 'YDir', 'normal');       % make (0,0) bottom-left like in your plot
set(gca, 'FontSize', 16);
box on

saveas(gcf, 'extra1_wswitch', 'png');


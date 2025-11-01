clear; clc; close all;

% ---- infinite square well energies (n = 1..5) ----
%E_num = [4.934802200545, 19.739208802179, 44.413219804902, 78.956835255528, 123.370055130650];   % numerical
%E_ana = [4.934802200545, 19.739208802179, 44.413219804902, 78.956835208715, 123.370055013617];   % analytical


% ---- parabolic well energies (n = 1..5) ----
E_num = [0.50, 1.50, 2.50, 3.50, 4.50];   % numerical
E_ana = [0.50, 1.50, 2.50, 3.50, 4.50];   % analytical

% ---- layout ----
xL = 0; xR = 1;                               % walls of the box
Eall = [E_num(:); E_ana(:)];
yr  = max(range(Eall), 1);                    % avoid zero span
ylo = min(Eall) - 0.22*yr;  yhi = max(Eall) + 0.22*yr;

dx_head = 0.035*(xR-xL);  dy_head = 0.06*(yhi-ylo);
shortFrac = 1.0;                              % length of BOTH ticks
xmid = 0.5*(xL+xR); halfW = 0.5*(xR-xL)*shortFrac;

figure(1); clf; hold on;

% walls + base + arrowheads
plot([xL xL],[ylo yhi],'k','LineWidth',3);
plot([xR xR],[ylo yhi],'k','LineWidth',3);
plot([xL xR],[ylo ylo],'k','LineWidth',3);          % horizontal x-axis line
patch([xL-dx_head xL xL+dx_head],[yhi-dy_head yhi yhi-dy_head],'k','EdgeColor','none');
patch([xR-dx_head xR xR+dx_head],[yhi-dy_head yhi yhi-dy_head],'k','EdgeColor','none');
text(xL,  yhi+0.04*yr, 'V = \infty', 'HorizontalAlignment','center');
text(xR,  yhi+0.04*yr, 'V = \infty', 'HorizontalAlignment','center');
%text(0.5, yhi+0.07*yr, 'Energy', 'HorizontalAlignment','center', 'FontWeight','bold');

% legend proxies (match styles)
hNum = plot(nan,nan,'b--','LineWidth',3);
hAna = plot(nan,nan,'k--','LineWidth',3);

% draw levels n = 1..5
for n = 1:5
  ya = E_ana(n);  yn = E_num(n);

  % analytical (black dotted), same length as numerical
  plot([xmid-halfW, xmid+halfW], [ya ya], 'k--', 'LineWidth', 3);

  % numerical (blue solid)
  plot([xmid-halfW, xmid+halfW], [yn yn], 'b-',  'LineWidth', 3);

  % right-side index labels E1..E5 at the analytical level
  text(xR + 0.06, ya, sprintf('E_%d', n), 'Interpreter','tex', ...
       'HorizontalAlignment','left', 'FontSize', 12);
end

% cosmetics
xlim([-0.15 1.35]); ylim([ylo-0.05*yr, yhi+0.10*yr]); axis off;
%lgd = legend([hNum hAna], {'E_{numerical}','E_{analytical}'}, ...
%             'Interpreter','tex', 'Location','north', 'FontSize', 14);

set(gcf,'Color','w');

print -dpng -r300 energy_pw.png


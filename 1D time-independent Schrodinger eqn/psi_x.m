clear; clc;

% ---- Load stitched modes (change names if yours differ) ----
m1 = load('fsw_n1.dat');   % ground state
m2 = load('fsw_n3.dat');   % 1st excited
m3 = load('fsw_n5.dat');   % 2nd excited
m4 = load('fsw_n7.dat');   % 3rd excited
%m5 = load('fsw_n9.dat');   % 4th excited

% Use column 5 (= |psi|^2) as in your current script
x1 = m1(:,1); psi1 = m1(:,3); V1 = m1(:,2);
x2 = m2(:,1); psi2 = m2(:,3); V2 = m2(:,2);
x3 = m3(:,1); psi3 = m3(:,3); V3 = m3(:,2);
x4 = m4(:,1); psi4 = m4(:,3); V4 = m4(:,2);
%x5 = m5(:,1); psi5 = m5(:,3); V5 = m5(:,2);

% Use the widest x-grid among the files for the potential outline
xV = x4; V = V4;

% ---- Normalize probability densities to max = 1 (for clean stacking) ----
normmax = @(v) v./max(abs(v));
psi1 = normmax(psi1);
psi2 = normmax(psi2);
psi3 = normmax(psi3);
psi4 = normmax(psi4);
%psi5 = normmax(psi5);

% ---- Choose vertical positions (display offsets, not actual energies) ----
E1 = 1.0; E2 = 3.0; E3 = 5.0; E4 = 7.0; %E5 = 9.0;

amp = 0.8;                      % visual amplitude scale for |psi|^2
y1 = E1 + amp*psi1;
y2 = E2 + amp*psi2;
y3 = E3 + amp*psi3;
y4 = E4 + amp*psi4;
%y5 = E5 + amp*psi5;

% ---- Scale the potential so it fits under the top state (purely cosmetic) ----
maxV      = max(V);
targetTop = E4 + 0.8;
s  = targetTop / maxV;
yV = s * V;

% ---- Plot ----
figure(1); clf; hold on;

% Potential outline (bold black)
plot(xV, yV, 'k-', 'LineWidth', 3.0, 'HandleVisibility','off');

% Energy guide lines (no legend entry)
plot([min(xV) max(xV)],[E1 E1],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E2 E2],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E3 E3],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E4 E4],'k:','LineWidth',1.2,'HandleVisibility','off');
%plot([min(xV) max(xV)],[E5 E5],'k:','LineWidth',1.2,'HandleVisibility','off');

% Wavefunctions (|psi|^2 stacks)
p1 = plot(x1, y1, 'r-', 'LineWidth', 2);   % Ground state
p2 = plot(x2, y2, 'b-', 'LineWidth', 2);   % 1st excited
p3 = plot(x3, y3, 'c-', 'LineWidth', 2);   % 2nd excited
p4 = plot(x4, y4, 'm-', 'LineWidth', 2);   % 3rd excited  (fixed from 'br-')
%p5 = plot(x5, y5, 'g-', 'LineWidth', 2);   % 4th excited

% ---- Axes styling: boxed and ticks inward ----
xlabel('x','FontSize',20,'FontWeight','bold');
ylabel('\psi(x)','FontSize',20,'FontWeight','bold');

set(gca, ...
    'Box','on', ...             % draw the full box
    'LineWidth',1.4, ...
    'TickDir','in', ...         % ticks inside the box
    'TickLength',[0.018 0.018], ...
    'XMinorTick','off','YMinorTick','off', ...
    'Layer','top', ...
    'FontSize',16);

% y-ticks as n-labels
yticks([E1 E2 E3 E4]);
yticklabels({'n=1','n=2','n=3','n=4'});

xlim([min(xV) max(xV)]);
ylim([0, targetTop + 0.4]);

% Optional legend
% lgd = legend([p1 p2 p3 p4 p5], ...
%     {'Ground State','1st Excited','2nd Excited','3rd Excited','4th Excited'}, ...
%     'Location','northeast','FontSize',14,'Box','off');

set(gcf,'color','w'); drawnow;

% ---- Save ----
print(gcf,'fsw_psi_boxed','-dpng','-r300');


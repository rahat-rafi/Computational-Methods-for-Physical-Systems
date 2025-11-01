clear; clc;

% ---- Load stitched modes (change names if yours differ) ----
m1 = load('fsw_n1.dat');   % ground state
m2 = load('fsw_n3.dat');   % 1st excited
m3 = load('fsw_n5.dat');   % 2nd excited
m4 = load('fsw_n7.dat');   % 3rd excited
m5 = load('fsw_n9.dat');   % 4th excited

x1 = m1(:,1); psi1 = m1(:,3); V1 = m1(:,2);
x2 = m2(:,1); psi2 = m2(:,3); V2 = m2(:,2);
x3 = m3(:,1); psi3 = m3(:,3); V3 = m3(:,2);
x4 = m4(:,1); psi4 = m4(:,3); V4 = m4(:,2);
x5 = m5(:,1); psi5 = m5(:,3); V5 = m5(:,2);

% Use the widest x-grid among the files for the potential outline
xV = x5; V = V5;

% ---- Normalize wavefunctions ----
psi1 = psi1 ./ max(abs(psi1));
psi2 = psi2 ./ max(abs(psi2));
psi3 = psi3 ./ max(abs(psi3));
psi4 = psi4 ./ max(abs(psi4));
psi5 = psi5 ./ max(abs(psi5));

% ---- Choose vertical positions (just display offsets, not actual energies) ----
E1 = 1.0;     % Ground shown at y = 1.0
E2 = 3.0;     % 1st excited shown at y = 2.5
E3 = 5.0;     % 2nd excited shown at y = 4.0
E4 = 7.0;     % Ground shown at y = 1.0
E5 = 9.0;     % 1st excited shown at y = 2.5

amp = 0.8;    % visual amplitude scale for wavefunctions
y1 = E1 + amp*psi1;
y2 = E2 + amp*psi2;
y3 = E3 + amp*psi3;
y4 = E4 + amp*psi4;
y5 = E5 + amp*psi5;

% ---- Scale the potential to fit nicely in the panel ----
maxV = max(V);                     % actual max of V over domain
targetTop = E5 + 0.8;              % put the box just under the top
s = targetTop / maxV;              % scale factor
yV = s * V;                        % scaled V(x) for display

% ---- Plot ----
figure(1); clf; hold on;

% Infinite square well outline (no legend entry)
plot(xV, yV, 'm--', 'LineWidth', 2, 'HandleVisibility', 'off');

% Energy guide lines (no legend entry)
plot([min(xV) max(xV)],[E1 E1],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E2 E2],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E3 E3],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E4 E4],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([min(xV) max(xV)],[E5 E5],'k:','LineWidth',1.2,'HandleVisibility','off');

% Wavefunctions (with legend)
p1 = plot(x1, y1, 'r-', 'LineWidth', 2);   % Ground state
p2 = plot(x2, y2, 'b-', 'LineWidth', 2);   % 1st excited
p3 = plot(x3, y3, 'k-', 'LineWidth', 2);   % 2nd excited
p4 = plot(x4, y4, 'y-', 'LineWidth', 2);   % 3rd state
p5 = plot(x5, y5, 'g-', 'LineWidth', 2);   % 4th excited

% Labels, ticks, styling
xlabel('L','FontSize',20,'FontWeight','bold');
ylabel('\psi(L)','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',16,'LineWidth',1.3);

% Replace y-axis ticks with n= labels at the three display levels
yticks([E1 E2 E3 E4 E5]);
yticklabels({'n=1','n=2','n=3', 'n=4','n=5'});

% Legend (only the states)
%legend([p1 p2 p3 p4 p5 p6 p7 p8], {'Ground State','1st Excited State','2nd Excited State', '3rd excited', '4', '5', '6', '7', '8'}, ...
%       'Location','northeast','FontSize',14);

xlim([min(xV) max(xV)]);
ylim([0, targetTop+0.4]);

set(gcf,'color','w');

% ---- Save figure as fig.png ----
print(gcf,'fsw_psi','-dpng','-r300');


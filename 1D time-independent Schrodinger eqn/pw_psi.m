clear; clc;

% ---- Load stitched modes (change names if yours differ) ----
m1 = load('pw_n0.dat');   % ground state
m2 = load('pw_n1.dat');   % 1st excited
m3 = load('pw_n2.dat');   % 2nd excited
m4 = load('pw_n3.dat');   % 3rd excited
m5 = load('pw_n4.dat');   % 4th excited

x1 = m1(:,1); psi1 = m1(:,3); V1 = m1(:,2);
x2 = m2(:,1); psi2 = m2(:,3); V2 = m2(:,2);
x3 = m3(:,1); psi3 = m3(:,3); V3 = m3(:,2);
x4 = m4(:,1); psi4 = m4(:,3); V4 = m4(:,2);
x5 = m5(:,1); psi5 = m5(:,3); V5 = m5(:,2);

% ---------------- crop to x in [xlo, xhi] ----------------
xlo = -5; xhi = 5;
mask1 = (x1>=xlo & x1<=xhi); x1 = x1(mask1); psi1 = psi1(mask1); V1 = V1(mask1);
mask2 = (x2>=xlo & x2<=xhi); x2 = x2(mask2); psi2 = psi2(mask2); V2 = V2(mask2);
mask3 = (x3>=xlo & x3<=xhi); x3 = x3(mask3); psi3 = psi3(mask3); V3 = V3(mask3);
mask4 = (x4>=xlo & x4<=xhi); x4 = x4(mask4); psi4 = psi4(mask4); V4 = V4(mask4);
mask5 = (x5>=xlo & x5<=xhi); x5 = x5(mask5); psi5 = psi5(mask5); V5 = V5(mask5);

% Use one grid for potential outline (after cropping)
xV = x5; V = V5;

% ---- Normalize wavefunctions *after* cropping ----
psi1 = psi1 ./ max(abs(psi1));
psi2 = psi2 ./ max(abs(psi2));
psi3 = psi3 ./ max(abs(psi3));
psi4 = psi4 ./ max(abs(psi4));
psi5 = psi5 ./ max(abs(psi5));

% ---- Display offsets (just for stacking) ----
E1 = 1.0; E2 = 3.0; E3 = 5.0; E4 = 7.0; E5 = 9.0;

amp = 0.8;                      % visual amplitude scale
y1 = E1 + amp*psi1;
y2 = E2 + amp*psi2;
y3 = E3 + amp*psi3;
y4 = E4 + amp*psi4;
y5 = E5 + amp*psi5;

% ---- Scale the potential to sit below the top state (use cropped V) ----
maxV = max(V);
targetTop = E5 + 2.0;
s = (targetTop) / maxV;
yV = s * V;

% ---- Plot ----
figure(1); clf; hold on;

% Potential outline
plot(xV, yV, 'k-', 'LineWidth', 3.5, 'HandleVisibility', 'off');

% Energy guide lines (no legend entry)
plot([xlo xhi],[E1 E1],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([xlo xhi],[E2 E2],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([xlo xhi],[E3 E3],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([xlo xhi],[E4 E4],'k:','LineWidth',1.2,'HandleVisibility','off');
plot([xlo xhi],[E5 E5],'k:','LineWidth',1.2,'HandleVisibility','off');

% Wavefunctions
p1 = plot(x1, y1, 'r-', 'LineWidth', 2);
p2 = plot(x2, y2, 'c-', 'LineWidth', 2);
p3 = plot(x3, y3, 'b-', 'LineWidth', 2);
p4 = plot(x4, y4, 'm-', 'LineWidth', 2);
p5 = plot(x5, y5, 'g-', 'LineWidth', 2);

% ---- Axes styling: boxed axes with ticks inward ----
xlabel('x','FontSize',20,'FontWeight','bold');
ylabel('\psi(x)','FontSize',20,'FontWeight','bold');

set(gca, ...
    'Box','on', ...             % draw full box
    'LineWidth',1.4, ...
    'TickDir','in', ...         % ticks inward
    'TickLength',[0.018 0.018], ...
    'XMinorTick','off','YMinorTick','off', ...
    'Layer','top', ...
    'FontSize',16);

% y-ticks as n-labels
yticks([E1 E2 E3 E4 E5]);
yticklabels({'n=1','n=2','n=3','n=4','n=5'});

xlim([xlo xhi]);
ylim([0, targetTop]);

set(gcf,'color','w'); drawnow;

% Save
print(gcf,'pw_psi_boxed','-dpng','-r300');


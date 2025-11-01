clear; clc;

% ---- Load stitched modes ----
m1 = load('isw_n1.dat');   % ground state
m2 = load('isw_n2.dat');   % 1st excited
m3 = load('isw_n3.dat');   % 2nd excited
m4 = load('isw_n4.dat');   % 3rd excited
m5 = load('isw_n5.dat');   % 4th excited

x1 = m1(:,1);  psi1 = m1(:,5);
x2 = m2(:,1);  psi2 = m2(:,5);
x3 = m3(:,1);  psi3 = m3(:,5);
x4 = m4(:,1);  psi4 = m4(:,5);
x5 = m5(:,1);  psi5 = m5(:,5);

% (Optional) restrict to the well domain [0,1] if files are wider
keep = @(x) (x>=0 & x<=1);
x1=x1(keep(x1)); psi1=psi1(keep(x1));
x2=x2(keep(x2)); psi2=psi2(keep(x2));
x3=x3(keep(x3)); psi3=psi3(keep(x3));
x4=x4(keep(x4)); psi4=psi4(keep(x4));
x5=x5(keep(x5)); psi5=psi5(keep(x5));

% ---- Normalize each state to max|psi|=1 (for clean stacking) ----
normmax = @(v) v./max(abs(v));
psi1 = normmax(psi1);
psi2 = normmax(psi2);
psi3 = normmax(psi3);
psi4 = normmax(psi4);
psi5 = normmax(psi5);

% ---- Display offsets (only for visualization) ----
E1 = 1.0; E2 = 2.0; E3 = 3.0; E4 = 4.0; E5 = 5.0;
amp = 0.38;                 % visual amplitude scale (fits nicely in the box)
y1 = E1 + amp*psi1;
y2 = E2 + amp*psi2;
y3 = E3 + amp*psi3;
y4 = E4 + amp*psi4;
y5 = E5 + amp*psi5;

% ---- Plot ----
figure(1); clf; hold on;

% Infinite-well walls (thick black at x=0 and x=1)
yl = [0, 5.9];                    % pad a bit above E5
plot([0 0], yl, 'k-', 'LineWidth', 2.8, 'HandleVisibility','off');
plot([1 1], yl, 'k-', 'LineWidth', 2.8, 'HandleVisibility','off');

% Energy guide lines (dotted; no legend)
plot([0 1],[E1 E1],'k:','LineWidth',1.1,'HandleVisibility','off');
plot([0 1],[E2 E2],'k:','LineWidth',1.1,'HandleVisibility','off');
plot([0 1],[E3 E3],'k:','LineWidth',1.1,'HandleVisibility','off');
plot([0 1],[E4 E4],'k:','LineWidth',1.1,'HandleVisibility','off');
plot([0 1],[E5 E5],'k:','LineWidth',1.1,'HandleVisibility','off');

% Wavefunctions
p1 = plot(x1, y1, 'r-', 'LineWidth', 2.0);
p2 = plot(x2, y2, 'b-', 'LineWidth', 2.0);
p3 = plot(x3, y3, 'c-', 'LineWidth', 2.0);
p4 = plot(x4, y4, 'm-', 'LineWidth', 2.0);
p5 = plot(x5, y5, 'g-', 'LineWidth', 2.0);

% ---- Axes styling: boxed and ticks inward (like the reference figure) ----
xlabel('x','FontSize',18,'FontWeight','bold');
ylabel('|\psi(x)|^2','FontSize',18,'FontWeight','bold');

set(gca, ...
    'Box','on', ...             % draw the full box
    'LineWidth',1.4, ...
    'TickDir','in', ...         % ticks inside the box
    'TickLength',[0.018 0.018], ...
    'XMinorTick','off', 'YMinorTick','off', ...
    'Layer','top', ...          % ticks/grid on top of plots
    'FontSize',14);

xlim([0 1]);
ylim(yl);

% Place y-ticks at the displayed "levels" with n-labels
yticks([E1 E2 E3 E4 E5]);
yticklabels({'n=1','n=2','n=3','n=4','n=5'});

% Optional: light vertical reference ticks (like in the sample)
xticks(0:0.2:1);

set(gcf,'color','w');
drawnow;
print(gcf,'isw_psi2_boxed','-dpng','-r300');


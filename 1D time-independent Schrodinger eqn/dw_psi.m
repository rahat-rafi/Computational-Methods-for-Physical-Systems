% plot_delta_vertical_spike.m
% Read dw_n0.dat and plot ψ(x) stacked at n=1.
% Potential is shown as a baseline with a single vertical downward spike at x=0.

clear; clc; close all;

% -------- Load & sort --------
D = load('dw_n0.dat');                 % cols: x V psi_num psi_exact ...
x          = D(:,1);
psi_num    = D(:,3);
psi_exact  = (size(D,2)>=4) * D(:,4);

[ x, ord ] = sort(x);
psi_num    = psi_num(ord);
if size(D,2)>=4, psi_exact = psi_exact(ord); else, psi_exact = []; end

% -------- Window & normalize --------
xwin = [-2, 2];
idx  = (x>=xwin(1)) & (x<=xwin(2));
x         = x(idx);
psi_num   = psi_num(idx);
if ~isempty(psi_exact), psi_exact = psi_exact(idx); end

psi_num = psi_num ./ max(abs(psi_num));
if ~isempty(psi_exact), psi_exact = psi_exact ./ max(abs(psi_exact)); end

% -------- Stacking level (n = 1) --------
E1  = 1.0;              % where we display the bound state
amp = 0.8;              % visual amplitude
y1  = E1 + amp*psi_num;
if ~isempty(psi_exact), y1e = E1 + amp*psi_exact; end

% -------- Delta potential as baseline + single vertical spike --------
baseY      = E1 - 0.10;      % baseline slightly below level line
spikeDepth = 0.55;           % how far the spike goes down from the baseline
x0         = 0.0;            % spike location

% -------- Plot --------
figure(1); clf; hold on;

% (1) Baseline
plot([xwin(1) xwin(2)], [baseY baseY], 'k-', 'LineWidth', 2.5, 'HandleVisibility','off');

% (2) Vertical downward spike at x=0
plot([x0 x0], [baseY baseY - spikeDepth], 'k-', 'LineWidth', 4, 'HandleVisibility','off');

% Optional tiny arrow head (comment out if you don't want it)
ah = 0.08;                     % arrowhead size
plot([x0-0.05 x0 x0+0.05], [baseY-spikeDepth+ah, baseY-spikeDepth, baseY-spikeDepth+ah], ...
     'k-', 'LineWidth', 3, 'HandleVisibility','off');

% (3) Energy guide (n=1)
plot([xwin(1) xwin(2)], [E1 E1], 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');

% (4) ψ overlays
hNum = plot(x, y1,  'r-',  'LineWidth', 2.4, 'DisplayName','\psi_{\rm numerical}(x)');
if ~isempty(psi_exact)
  hAna = plot(x, y1e, 'b--', 'LineWidth', 2.0, 'DisplayName','\psi_{\rm exact}(x)');
end

% -------- Axes & labels --------
xlabel('x','FontSize',20,'FontWeight','bold');
ylabel('\psi(x)','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',16,'LineWidth',1.3,'TickDir','in','Box','on');
yticks([E1]); yticklabels({'n=1'});
xlim(xwin);
ylim([0, E1+1.2]);


% Legend
%if exist('hAna','var')
%  legend([hNum hAna], 'Location','northeast','Box','off','FontSize',14);
%else
%  legend(hNum, 'Location','northeast','Box','off','FontSize',14);
%end

set(gcf,'color','w');
print(gcf,'dw_psi','-dpng','-r300');


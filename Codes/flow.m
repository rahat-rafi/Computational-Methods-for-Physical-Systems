clear all;
data = load('case1_field_200K_step10.dat');
quiver(data(:,1),data(:,2),data(:,3),data(:,4),'b')
hold on
rectangle('position', [401 1 50 50], 'facecolor', 'r','edgecolor', 'r')
axis([0 870 1 451])
set(gca,'FontSize',16)
hold on
saveas(gcf,'Bc1_200','epsc')

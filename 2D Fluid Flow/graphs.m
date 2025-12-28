clear all
x1 = load('case2_psi_ix390_snap.dat');

n = 449;

figure;
for i = 1:2:10
   plot(x1((i-1)*n+1:i*n,2),x1((i-1)*n+1:i*n,6))
   hold on
end
hold off
xlabel('y','fontsize',2);
ylabel('\Omega', 'Fontsize',3);

hold off
axis([0 450 -0.013 0.001])
set(gca,'FontSize',16)
hold on
saveas(gcf,'Bc2_200_omegavsy','epsc')

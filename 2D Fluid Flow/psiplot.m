clear all
x1 = load('case1_field.dat');

n = 449;

figure;
for i = 1:2:10
   plot(x1((i-1)*n+1:i*n,2),x1((i-1)*n+1:i*n,6))
   hold on
end
hold off
xlabel('y','fontsize',2);
ylabel('u', 'Fontsize',3);

hold off
axis([0 450 0 1.2])
set(gca,'FontSize',16)
hold on
saveas(gcf,'Bc2_200_u','epsc')

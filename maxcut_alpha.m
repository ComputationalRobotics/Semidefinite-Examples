clc; clear; close all;

alpha = 0.878;

t = -1:0.01:1;

lb = alpha * (1-t);
ub = 2/pi * acos(t);

figure;
plot(t,ub,'LineWidth',4);
hold on
plot(t,lb,'LineWidth',4);
axis equal
grid on
ax = gca;
ax.FontSize = 16;
legend('$2/\pi \arccos t$', '$\alpha (1 - t)$','FontSize',20,'Interpreter','latex');
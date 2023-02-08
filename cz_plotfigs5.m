clc
clear all
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

% Load Data
load('./C_EPMC_nh33_P5.mat', 'qts', 'qdts', 'UxwC');

largs = {'-', 'LineWidth', 2};

[~, ki] = max(UxwC(end-1,:));  % point of max frequency
% ki = find(UxwC(end-2,:)./(2*UxwC(end-1,:))<1, 1 );  % find point when this mode becomes underdamped

figure(1)
clf()
subplot(2,1,1)
semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
ylabel('Frequency (rad/s)')
ax = axes('Position', [0.575 0.65 0.3 0.15]);
semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
xlim(10.^UxwC(end,[ki-1 end]))

subplot(2,1,2)
semilogx(10.^UxwC(end,:), UxwC(end-2,:)./(2*UxwC(end-1,:))*100, largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-2,ki-1)./(2*UxwC(end-1,ki-1))*100, 'ro', 'MarkerFaceColor', 'r'); hold on
ylabel('Eff. Damping (\%)')
xlabel('Modal Amplitude $q_1$')
ax = axes('Position', [0.575 0.25 0.3 0.15]);
semilogx(10.^UxwC(end,:), UxwC(end-2,:)./(2*UxwC(end-1,:))*100, largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-2,ki-1)./(2*UxwC(end-1,ki-1))*100, 'ro', 'MarkerFaceColor', 'r'); hold on
xlim(10.^UxwC(end,[ki-1 end]))
set(gcf, 'Color', 'white')

figure(2)
clf()
surf(qts([1:end 1],1:ki,1), qts([1:end 1],1:ki,2), qdts([1:end 1],1:ki,1), 'EdgeColor', 'none'); hold on
plot3(qts([1:end 1],1:4:ki,1), qts([1:end 1],1:4:ki,2), qdts([1:end 1],1:4:ki,1), 'k-')
colormap(jet)
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
set(gca, 'View', [15 5])

figure(3)
clf()
surf(qts([1:end 1],ki:end,1), qts([1:end 1],ki:end,2), qdts([1:end 1],ki:end,1), 'EdgeColor', 'none'); hold on
plot3(qts([1:end 1],ki:4:end,1), qts([1:end 1],ki:4:end,2), qdts([1:end 1],ki:4:end,1), 'k-')
colormap(jet)
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
set(gca, 'View', [-88 10])
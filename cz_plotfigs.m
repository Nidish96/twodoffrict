clc
clear all
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

% Load Data
load('./C_EPMC_nh33_P4.mat', 'qts', 'qdts', 'UxwC');

largs = {'-', 'LineWidth', 2};

[~, ki] = max(UxwC(end-2,:));  % Plot until point of max dissipation

figure(1)
clf()
subplot(2,1,1)
semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
ylabel('Frequency (rad/s)')

subplot(2,1,2)
semilogx(10.^UxwC(end,:), UxwC(end-2,:)./(2*UxwC(end-1,:))*100, largs{:}); hold on
semilogx(10.^UxwC(end,ki-1), UxwC(end-2,ki-1)./(2*UxwC(end-1,ki-1))*100, 'ro', 'MarkerFaceColor', 'r'); hold on

ylabel('Eff. Damping (\%)')
xlabel('Modal Amplitude $q_1$')
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
clc
clear all
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

% Load Data
load('./DATS/C_EPMC_nh33_P5.mat', 'uts', 'udts', 'qts', 'qdts', 'UxwC', 'Fts');
h = 1:2:33;
Nt = 1024;
t = linspace(0, 2*pi, Nt+1);  t = t(1:Nt)';

Nhc = sum((h==0)+2*(h~=0));

%% Compute Invariant Manifold
Phis = [1 1;1 -1]/sqrt(2);  % Mode shapes used by Prof. Quinn
U1 = UxwC(1:2:end-3,:).*10.^UxwC(end,:);                           U2 = UxwC(2:2:end-3,:).*10.^UxwC(end,:);
Q1 = kron(eye(Nhc), Phis(:,1)')*UxwC(1:end-3,:).*10.^UxwC(end,:);  Q2 = kron(eye(Nhc), Phis(:,2)')*UxwC(1:end-3,:).*10.^UxwC(end,:);

uts = zeros(Nt, size(UxwC,2)+1, 2);
uts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) U1], 0);
uts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) U2], 0);

udts = zeros(Nt, size(UxwC,2)+1, 2);
udts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) U1], 1);
udts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) U2], 1);

qts = zeros(Nt, size(UxwC,2)+1, 2);
qts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Q1], 0);
qts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Q2], 0);

qdts = zeros(Nt, size(UxwC,2)+1, 2);
qdts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Q1], 1);
qdts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Q2], 1);

%%
largs = {'-', 'LineWidth', 2};

[~, ki] = max(UxwC(end-1,:));  % point of max frequency
% ki = find(UxwC(end-2,:)./(2*UxwC(end-1,:))<1, 1 );  % find point when this mode becomes underdamped
ki = ki+20;

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
surf(qts([1:end 1],1:ki,1), qdts([1:end 1],1:ki,1), qts([1:end 1],1:ki,2), 'EdgeColor', 'none'); hold on
plot3(qts([1:end 1],1:4:ki,1), qdts([1:end 1],1:4:ki,1), qts([1:end 1],1:4:ki,2), 'k-')
colormap(jet)
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
set(gca, 'View', [15 5])

% kmax = size(UxwC,2); % MODIFY THIS TO CONTROL THE MAX AMPLITUDE TO PLOT
kmax = ki + 10;
figure(3)
clf()
surf(qts([1:end 1],ki:kmax,1), qdts([1:end 1],ki:kmax,1), qts([1:end 1],ki:kmax,2), 'EdgeColor', 'none'); hold on
plot3(qts([1:end 1],ki:4:kmax,1), qdts([1:end 1],ki:4:kmax,1), qts([1:end 1],ki:4:kmax,2), 'k-')
colormap(jet)
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
set(gca, 'View', [-88 10])

figure(4)
clf()
subplot(3,1,1)
plot(t, uts(:, ki, 2), 'LineWidth', 2)
grid on
xlabel('Scaled Time')
ylabel('Disp. $x_2$')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi,'XTickLabel',{'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'},'TickLabelInterpreter','latex')
subplot(3,1,2)
plot(t, udts(:, ki, 2), 'LineWidth', 2)
grid on
xlabel('Scaled Time')
ylabel('Velocity $\dot{x}_2$')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi,'XTickLabel',{'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'},'TickLabelInterpreter','latex')
subplot(3,1,3)
plot(t, Fts(:, ki, 2), 'LineWidth', 2)
grid on
xlabel('Scaled Time')
ylabel('Force $F_2$')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi,'XTickLabel',{'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'},'TickLabelInterpreter','latex')
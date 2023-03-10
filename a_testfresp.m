% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfig = true;
%% Parameters
bt = 0.1;
nlpars = [1;0;0.5];  % [alpha; delta; muN]
p = struct('M', eye(2), 'C', bt*[0 0;0 1], 'K', [2 -1;-1 2], ...
    'nlpars', nlpars, 'fv', [1;0]*10, 'kt', inf);  % set kt to inf for rigid coulomb
epN = 1e3;

%% Transient Simulation Using Implicit Acceleration (See NEWMARKBMARCH for modifying parameters)
Om = 0.9;
Ncyc = 50;
T = 2*pi/Om*Ncyc;
fsamp = 10;

t = (0:1/fsamp:T);
[t, ut, udt, uddt, fnlt] = NEWMARKBMARCH(p.fv.*cos(Om*t), zeros(2,1), zeros(2,1), t, p, epN);
ut = ut'; udt = udt'; uddt = uddt'; fnlt = fnlt';
%% HB Simulation
h = 1:2:33;  % Choose harmonics to balance. Ex: h = [0:5];
[~,~, hinds,rinds,iinds] = HINDS(1, h);
Nhc = sum((h==0)+2*(h~=0));
Nt = 1024;

Fl = zeros(2*Nhc, 1);
Fl(rinds(1:2)) = p.fv;
U0 = HARMONICSTIFFNESS(p.M, p.C, p.K, Om, h)\Fl;
% U0(2:2:end) = 0;

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
UH = fsolve(@(U) HBRESFUN([U;Om], Fl, h, Nt, p, epN), U0, opt);
[~, ~, ~, FNLH] = HBRESFUN([UH;Om], Fl, h, Nt, p, epN);

D1 = HARMONICSTIFFNESS(0, 1, 0, Om, h);
th = linspace(0, 2*pi/Om, Nt+1); th = th(1:end-1);
uh = TIMESERIES_DERIV(Nt, h, reshape(UH, 2, Nhc)', 0);
udh= TIMESERIES_DERIV(Nt, h, D1*reshape(UH, 2, Nhc)', 0);
fnlh=TIMESERIES_DERIV(Nt, h, reshape(FNLH, 2, Nhc)', 0);

%% Plot Results
figure(1)
clf()
subplot(3,2,1)
plot(ut(:,1), ut(:,2), 'k-'); hold on
plot(uh(:,1), uh(:,2), 'LineWidth', 2)
xlabel('Disp $x_1$ (m)')
ylabel('Disp $x_2$ (m)')
title('\textbf{Conf. Space}')

subplot(3,2,3)
UH1 = UH(1:2:end);
stem(nan, nan); hold on
stem(h*Om, abs([UH1(hinds); UH1(rinds)-1j*UH1(iinds)]), '-', 'filled')

xlabel('Frequency (rad/s)')
ylabel('Disp $x_1$ (m)')

subplot(3,2,2)
plot(ut(:,1), fnlt(:,1), 'k.'); hold on
plot(uh(:,1), fnlh(:,1), '.')
title('\textbf{Transient Results}')
xlabel('Displacement $x_1$ (m)')
ylabel('NL force 1 (N)')

subplot(3,2,4)
if isinf(p.kt)
    plot(udt(:,2), fnlt(:,2), 'k.'); hold on
    plot(udh(:,2), fnlh(:,2), '.')
    legend('Tr', 'HB', 'Location', 'northwest')
    xlabel('Velocity $v_2$ (m/s)')
else
    plot(ut(:,2), fnlt(:,2), 'k.'); hold on
    plot(uh(:,2), fnlh(:,2), '.')
    legend('Tr', 'HB', 'Location', 'northwest')
    xlabel('Displacement $x_2$ (m)')
end
ylabel('NL force 2 (N)')

subplot(3,2,5)
plot(kron((0:Ncyc-1)*2*pi/Om,ones(1,Nt))+repmat(th, 1, Ncyc), repmat(uh(:,1), Ncyc, 1), 'k.-'); hold on
plot(t, ut(:,1), 'LineWidth', 2); hold on
xlabel('Time (s)')
ylabel('Disp $x_1$ (m)')
xlim(T(end)-2*pi/Om*[4 0])

subplot(3,2,6)
plot(kron((0:Ncyc-1)*2*pi/Om,ones(1,Nt))+repmat(th, 1, Ncyc), repmat(uh(:,2), Ncyc, 1), 'k.-'); hold on
plot(t, ut(:,2), 'LineWidth', 2); hold on
xlabel('Time (s)')
ylabel('Disp $x_2$ (m)')
xlim(T(end)-2*pi/Om*[4 0])
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/A_TrComp_%d.png', h(end)), '-dpng')
end

%%
figure(2)
clf()
plot3(ut(:,1), ut(:,2), udt(:,1), 'k-'); hold on
plot3(uh(:,1), uh(:,2), udh(:,1), 'LineWidth', 4)
box on
legend('Transient', 'HBM', 'Location', 'best')
xlabel('Disp $x_1$ (m)')
ylabel('Disp $x_2$ (m)')
zlabel('Velocity $v_1$ (m)')
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/A_TrCompSS_%d.png', h(end)), '-dpng')
end

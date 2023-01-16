% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%% Parameters
bt = 0.1;
nlpars = [0.005;1e-10];  % [alpha, mu]. [1;2]
p = struct('M', eye(2), 'C', bt*[0 0;0 1], 'K', [2 -1;-1 2], ...
    'nlpars', nlpars, 'fv', [1;0]*10, 'epN', 1e1);

%% Linear Modes
[V, D] = eig(p.K, p.M);
[Wsr, si] = sort(sqrt(diag(D)));
V = V(:, si);

%% EPMC Setup
mi = 1;  % Mode of interest

As = -4;
Ae = -1;
da = 0.1;

h = [1:2:5];
Nt = 1024;

Nhc = sum((h==0)+2*(h~=0));
Fl = zeros(2*Nhc, 1);
if h(1)==0
    Fl(3:4) = [1;0];
    U0 = [zeros(2,1); V(:,mi); V(:,mi); zeros(2*(Nhc-3),1)];
else
    Fl(1:2) = [1;0];
    U0 = [V(:,mi); V(:,mi); zeros(2*(Nhc-2),1)];
end

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
Uxws = fsolve(@(Uxw) EPMCRESFUN([Uxw; As], Fl, h, Nt, p), [U0;0;Wsr(mi)], opt);

%%
As = -20;
Ae = -1;
da = 0.1;
U0 = Uxws;
Sopt = struct('Nmax', 100, 'parametrization', 'orthogonal', 'dynamicDscale', 1);
UxwC = solve_and_continue(Uxws, @(Uxwa) EPMCRESFUN(Uxwa,Fl,h,Nt,p), As, Ae, da, Sopt);

figure(1)
clf()
subplot(2,1,1)
semilogx(10.^UxwC(end,:), UxwC(end-1,:), '.-')
ylabel('Frequency (rad/s)')
subplot(2,1,2)
semilogx(10.^UxwC(end,:), UxwC(end-2,:), '.-')
ylabel('xi')
xlabel('Amplitude')

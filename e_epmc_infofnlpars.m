% This conducts EPMC and plots the influence of the two parameters in the
% problem

% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfig = true;
%% Parameters
bt = 0.1;
nlpars = [1;0;0.5];  % [alpha, mu]
p = struct('M', eye(2), 'C', bt*[0 0;0 1], 'K', [2 -1;-1 2], ...
    'nlpars', nlpars, 'fv', [1;0]*10, 'kt', inf);
epN = 1e0;

%% Linear Modes
[V, Wsr] = eig(p.K, p.M);
[Wsr, si] = sort(sqrt(diag(Wsr)));
V = V(:, si);

% Fully Stuck Configuration
Lst = [1;0];
[Vst, Wst] = eig(Lst'*p.K*Lst, Lst'*p.M*Lst);
[Wst, si] = sort(sqrt(diag(Wst)));
Vst = Lst*Vst(:,si);

%% Setup 
h = 1:2:13;
Nt = 512;
Nhc = sum((h==0)+2*(h~=0));
[~,~,hinds,rinds,iinds] = HINDS(2,h);
Fl = zeros(2*Nhc, 1);
Fl(rinds(1)) = 1;

% Fully stuck initial guess
U0 = zeros(2*Nhc, 1);
U0([rinds(1:2) iinds(1:2)]) = repmat(Vst,2,1)/sqrt(2);
Uxw0 = [U0;Vst'*p.C*Vst;Wst];

% Influence of mu N
PARS = [1 0.25;
    1 0.5;
    1 0.75];
suff = 'muvar';
% % Influence of alpha
% PARS = [0.5 0.5;
%     1 0.5;
%     1.5 0.5];
% suff = 'alvar';

UxwC = cell(size(PARS,1),1);
for i=1:size(PARS,1)
    %% EPMC Setup        
    opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, ...
        'Display', 'iter', 'CheckGradients', false);
    p.nlpars([1 3]) = PARS(i,:)';
    
    Uxws = fsolve(@(Uxw) EPMCRESFUN([Uxw; -1], Fl, h, Nt, p, epN), Uxw0, opt);
    
    %% Conduct Continuation
    As = -1;
    Ae = 1;
    da = 0.025;  % 0.025
    U0 = Uxws;
    
    Copt = struct('Nmax', 100, 'solverchoice', 1, 'DynDscale', 1);
    Copt.Dscale = [zeros(Nhc*2,1); 0.05; Wst; 1.0];
    Copt.Dscale(hinds) = 1e-6;
    Copt.Dscale([rinds iinds]) = 1e0;
    UxwC{i} = CONTINUE(@(Uxwa) EPMCRESFUN(Uxwa,Fl,h,Nt,p,epN), Uxws, As, Ae, da, Copt);
end

%% Plot BackBones
largs = {'-', 'LineWidth', 2};

figure(7)
clf()
aa = gobjects(size(PARS,1),1);
for i=1:size(PARS,1)
    subplot(2,1,1)
    aa(i) = semilogx(10.^UxwC{i}(end,:), UxwC{i}(end-1,:), largs{:}); hold on
    legend(aa(i), sprintf('$\\alpha$ = %.1f N $m^{-3}$; $\\mu N$ = %.2f N', PARS(i,1), PARS(i,2)), 'interpreter', 'latex');

    subplot(2,1,2)
    semilogx(10.^UxwC{i}(end,:), UxwC{i}(end-2,:)./(2*UxwC{i}(end-1,:))*100, largs{:}); hold on
end
subplot(2,1,1)
legend(aa, 'Location', 'northwest', 'NumColumns', 1)
ylabel('Frequency (rad/s)')
subplot(2,1,2)
ylabel('Eff. Damping (\%)')
xlabel('Modal Amplitude $q_1$')
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/E_EPMCBBPCOMP_%d_%s.png', h(end), suff), '-dpng')
end
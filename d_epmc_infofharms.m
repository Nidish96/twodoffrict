% This conducts EPMC and plots the influence of different harmonic
% truncations

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
Hs = 1:4:33;
Nt = 512;
UxwC = cell(size(Hs));
for hi=1:length(Hs)
    %% EPMC Setup
    h = 1:2:Hs(hi);
    
    Nhc = sum((h==0)+2*(h~=0));
    [~,~,hinds,rinds,iinds] = HINDS(2,h);
    Fl = zeros(2*Nhc, 1);
    Fl(rinds(1)) = 1;
    
    % Fully stuck initial guess
    U0 = zeros(2*Nhc, 1);
    U0([rinds(1:2) iinds(1:2)]) = repmat(Vst,2,1)/sqrt(2);
    Uxw0 = [U0;Vst'*p.C*Vst;Wst];
    
    opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, ...
        'Display', 'iter', 'CheckGradients', false);
    
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
    UxwC{hi} = CONTINUE(@(Uxwa) EPMCRESFUN(Uxwa,Fl,h,Nt,p,epN), Uxws, As, Ae, da, Copt);
end

%% Plot BackBones
largs = {'-', 'LineWidth', 2};
figure(6)
clf()
aa = gobjects(size(Hs));
for hi=1:length(Hs)
    subplot(2,1,1)
    semilogx(10.^UxwC{hi}(end,:), UxwC{hi}(end-1,:), largs{:}); hold on

    subplot(2,1,2)
    aa(hi) = semilogx(10.^UxwC{hi}(end,:), UxwC{hi}(end-2,:)./(2*UxwC{hi}(end-1,:))*100, largs{:}); hold on
    legend(aa(hi), sprintf('H = %d', Hs(hi)));
end
subplot(2,1,1)
ylabel('Frequency (rad/s)')
ax = axes('Position', [0.2 0.75 0.4 0.15]);
for hi=1:length(Hs)
    semilogx(10.^UxwC{hi}(end,:), UxwC{hi}(end-1,:), largs{:}); hold on
    xlim([0.54 0.8])
end
subplot(2,1,2)
legend(aa, 'Location', 'northeast', 'NumColumns', 2, 'FontSize', 9)
ylabel('Eff. Damping (\%)')
xlabel('Modal Amplitude $q_1$')
set(gcf, 'Color', 'white')
if savfig
    export_fig('./FIGS/D_EPMCBBHCOMP.png', '-dpng')
end
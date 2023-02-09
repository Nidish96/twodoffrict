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
% [alpha, gamma, mu]
nlpars = [1 0 0.5;
    0 0 0.5;
    1 0 0;
    0 0 0.1;
    0 0 0.1];
nlpi = 5;
nlpars = nlpars(nlpi,:)';
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

%% EPMC Setup
h = 1:2:33;
Nt = 1024;

Nhc = sum((h==0)+2*(h~=0));
[~,~,hinds,rinds,iinds] = HINDS(2,h);
Fl = zeros(2*Nhc, 1);
Fl(rinds(1)) = 1;

% % Fully slipped initial guess
% mi = 1;  % Mode of interest
% U0 = zeros(2*Nhc, 1);
% U0([rinds(1:2) iinds(1:2)]) = repmat(V(:,mi),2,1)/sqrt(2);
% Uxw0 = [U0;V(:,mi)'*p.C*V(:,mi);Wsr(mi)];

% Fully stuck initial guess
U0 = zeros(2*Nhc, 1);
U0([rinds(1:2) iinds(1:2)]) = repmat(Vst,2,1)/sqrt(2);
Uxw0 = [U0;Vst'*p.C*Vst;Wst];

if nlpi==5  % Slipped initial guess
    U0([rinds(1:2) iinds(1:2)]) = repmat(V(:,1),2,1)/sqrt(2);
    Uxw0 = [U0;V(:,1)'*p.C*V(:,1);Wsr(1)];
end

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, ...
    'Display', 'iter', 'CheckGradients', false);
As = -2;
Ae = 1;
if nlpi==5
    Uxws = fsolve(@(Uxw) EPMCRESFUN([Uxw; Ae], Fl, h, Nt, p, epN), Uxw0, opt);
else
    Uxws = fsolve(@(Uxw) EPMCRESFUN([Uxw; As], Fl, h, Nt, p, epN), Uxw0, opt);
end

Ctran = zeros(length(h), Nhc);
[ci0,cin,hi0,hic,his] = HINDS(1,h);
Ctran(ci0,hi0) = eye(1);
Ctran(cin,hic) = eye(length(cin));
Ctran(cin,his) = 1j*eye(length(cin));

%% Conduct Continuation
As = -1;
Ae = 1;
if nlpi==2
    da = 0.0125;
elseif nlpi==4
    As = -2;
    da = 0.005;
elseif nlpi==5
    As = -4;
    da = 0.025;    
else
    da = 0.025;  % 0.025
end
U0 = Uxws;

Copt = struct('Nmax', 200, 'solverchoice', 1, 'DynDscale', 1);
if any(nlpi==[4 5])
    Copt.solverchoice = 2;
    Copt.Nmax = 400;
end
Copt.Dscale = [zeros(Nhc*2,1); 0.05; Wst; 1.0];
Copt.Dscale(hinds) = 1e-6;
Copt.Dscale([rinds iinds]) = 1e0;
if nlpi==5
    UxwC = CONTINUE(@(Uxwa) EPMCRESFUN(Uxwa,Fl,h,Nt,p,epN), Uxws, Ae, As, da, Copt);
    UxwC = UxwC(:, end:-1:1);
else
    UxwC = CONTINUE(@(Uxwa) EPMCRESFUN(Uxwa,Fl,h,Nt,p,epN), Uxws, As, Ae, da, Copt);
end

%% Save Force Harmonics
Fnls = zeros(Nhc*2, size(UxwC,2));
Fts = zeros(Nt, size(UxwC,2)+1, 2);
for ia=1:size(UxwC,2)
    [~, ~, ~, Fnls(:,ia)] = EPMCRESFUN(UxwC(:,ia),Fl,h,Nt,p,epN);
end
Fts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Fnls(1:2:end,:)], 0);
Fts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Fnls(2:2:end,:)], 0);

%% Invariant Manifold
uts = zeros(Nt, size(UxwC,2)+1, 2);
uts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) UxwC(1:2:end-3,:)], 0).*[0 10.^UxwC(end-1,:)];
uts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) UxwC(2:2:end-3,:)], 0).*[0 10.^UxwC(end-1,:)];

udts = zeros(Nt, size(UxwC,2)+1, 2);
udts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) UxwC(1:2:end-3,:)], 1).*[0 UxwC(end-1,:).*(10.^UxwC(end-1,:))];
udts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) UxwC(2:2:end-3,:)], 1).*[0 UxwC(end-1,:).*(10.^UxwC(end-1,:))];

Phis = [1 1;1 -1]/sqrt(2);  % Mode shapes used by Prof. Quinn
qts = zeros(Nt, size(UxwC,2)+1, 2);
qts(:,:,1) = Phis(1,1)*uts(:,:,1)+Phis(2,1)*uts(:,:,2);
qts(:,:,2) = Phis(1,2)*uts(:,:,1)+Phis(2,2)*uts(:,:,2);

qdts = zeros(Nt, size(UxwC,2)+1, 2);
qdts(:,:,1) = Phis(1,1)*udts(:,:,1)+Phis(2,1)*udts(:,:,2);
qdts(:,:,2) = Phis(1,2)*udts(:,:,1)+Phis(2,2)*udts(:,:,2);

%% 
largs = {'-', 'LineWidth', 2};
ki = min(50, size(UxwC,2));  % Maximum number of points to plot (see figure 4 below for reference)
if nlpi==4
    [~, ki] = max(UxwC(end-2,:));
elseif nlpi==5
    [~, ki] = max(UxwC(end-1,:));
end

figure(4)
clf()
subplot(2,1,1)
semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
if any(nlpi==[1 4 5])
    semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
end
ylabel('Frequency (rad/s)')
if nlpi==1
    ax = axes('Position', [0.2 0.75 0.4 0.15]);
    semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
    semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
    xlim(10.^UxwC(end,[max(1,ki-20) min(size(UxwC,2),ki+10)]))
elseif nlpi==5
    ax = axes('Position', [0.575 0.65 0.3 0.15]);
    semilogx(10.^UxwC(end,:), UxwC(end-1,:), largs{:}); hold on
    semilogx(10.^UxwC(end,ki-1), UxwC(end-1,ki-1), 'ro', 'MarkerFaceColor', 'r'); hold on
    xlim(10.^UxwC(end,[ki-1 end]))
end

subplot(2,1,2)
semilogx(10.^UxwC(end,:), UxwC(end-2,:)./(2*UxwC(end-1,:))*100, largs{:}); hold on
if any(nlpi==[1 4 5])
    semilogx(10.^UxwC(end,ki-1), UxwC(end-2,ki-1)./(2*UxwC(end-1,ki-1))*100, 'ro', 'MarkerFaceColor', 'r'); hold on
end
ylabel('Eff. Damping (\%)')
xlabel('Modal Amplitude $q_1$')
if nlpi==5
    ax = axes('Position', [0.575 0.25 0.3 0.15]);
    semilogx(10.^UxwC(end,:), UxwC(end-2,:)./(2*UxwC(end-1,:))*100, largs{:}); hold on
    semilogx(10.^UxwC(end,ki-1), UxwC(end-2,ki-1)./(2*UxwC(end-1,ki-1))*100, 'ro', 'MarkerFaceColor', 'r'); hold on
    xlim(10.^UxwC(end,[ki-1 end]))
end
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/C_EPMCBB_%d_P%d.png', h(end), nlpi), '-dpng');
end

figure(5)
clf()
surf(qts([1:end 1],1:ki,1), qts([1:end 1],1:ki,2), qdts([1:end 1],1:ki,1), 'EdgeColor', 'none'); hold on
if any(nlpi==[4 5])
    plot3(qts([1:end 1],1:4:ki,1), qts([1:end 1],1:4:ki,2), qdts([1:end 1],1:4:ki,1), 'k-')
else
    plot3(qts([1:end 1],1:ki,1), qts([1:end 1],1:ki,2), qdts([1:end 1],1:ki,1), 'k-')
end
colormap(jet)
switch nlpi
    case 1
        set(gca, 'View', [-30 15])
    case 2
        set(gca, 'View', [-40 15])
    case 5
        set(gca, 'View', [15 5])
end
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/C_EPMCIM_%d_P%d.png', h(end), nlpi), '-dpng');
end

figure(6)
clf()
surf(qts([1:end 1],ki:end,1), qts([1:end 1],ki:end,2), qdts([1:end 1],ki:end,1), 'EdgeColor', 'none'); hold on
if any(nlpi==[4 5])
    plot3(qts([1:end 1],ki:4:end,1), qts([1:end 1],ki:4:end,2), qdts([1:end 1],ki:4:end,1), 'k-')
else
    plot3(qts([1:end 1],ki:end,1), qts([1:end 1],ki:end,2), qdts([1:end 1],ki:end,1), 'k-')
end
colormap(jet)
switch nlpi
    case 1
        set(gca, 'View', [-30 15])
    case 2
        set(gca, 'View', [-40 15])
    case 5
        set(gca, 'View', [-88 10])
end
xlabel('Disp $q_1$');
ylabel('Velocity $d q_1/dt$')
zlabel('Disp $q_2$')
grid on
set(gcf, 'Color', 'white')
if savfig
    export_fig(sprintf('./FIGS/C_EPMCIM_%d_P%d_eh.png', h(end), nlpi), '-dpng');
end

save(sprintf('./DATS/C_EPMC_nh%d_P%d.mat', h(end), nlpi), 'uts', 'udts', 'qts', 'qdts', 'UxwC', 'Fnls', 'Fts');
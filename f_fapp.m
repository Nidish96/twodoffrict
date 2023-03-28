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

%% 
h = 1:2:33;
Nt = 1024;

Nhc = sum((h==0)+2*(h~=0));
[~,~,hinds,rinds,iinds] = HINDS(2,h);
Fl = zeros(2*Nhc, 1);
Fl(rinds(1:2)) = V(:,1);

U0w = zeros(2*Nhc+1,1);
U0w(iinds(1:2)) = V(:,1);
U0w(end) = Wsr(1);
% U0w([rinds(1:2) iinds(1:2)]) = repmat(Vst, 2,1);
% U0w(end) = Wst;

Fs = -4;
Fe = 3;

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
Uws = fsolve(@(Uw) FAPPRESFUN([Uw;Fs], Fl, h, Nt, p, epN), U0w, opt);

%%
Fs = -4;
Fe = 3;
df = 0.01;  % 0.05

Copt = struct('Nmax', 100, 'DynDscale', 1);
Copt.Dscale = [ones(2*Nhc,1)*max(abs(U0w(1:end-1))); 1.0; 1.0];

UwfC = CONTINUE(@(Uwf) FAPPRESFUN(Uwf, Fl, h, Nt, p, epN), Uws, Fs, Fe, df, Copt);

%% Calculate Force Harmonics
Fnls = zeros(Nhc*2, size(UwfC,2));
Fts = zeros(Nt, size(UwfC,2)+1, 2);
for ia=1:size(UwfC,2)
    [~, ~, ~, Fnls(:,ia)] = FAPPRESFUN(UwfC(:,ia),Fl,h,Nt,p,epN);
end
Fts(:,:,1) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Fnls(1:2:end,:)], 0);
Fts(:,:,2) = TIMESERIES_DERIV(Nt, h, [zeros(Nhc,1) Fnls(2:2:end,:)], 0);

%% Save Data
% save(sprintf('./DATS/F_FAPP_nh%d_P%d.mat', h(end), nlpi), 'UwfC', 'Fnls', 'Fts');

%% Plot Forces
figure(1)
clf()
% hold on
semilogx(10.^UwfC(end,:), UwfC(end-1,:), '.-')
xlabel('Force (N)')
ylabel('Natural Frequency (rad/s)')
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
epN = 1e0;

[V,Wsr] = eig(p.K, p.M);
[Wsr, si] = sort(sqrt(diag(Wsr)));
V = V(:,si);

Lst = [1;0];

%% HB Simulation
h = 1:2:33;  % Choose harmonics to balance. Ex: h = [0:5];
Nhc = sum((h==0)+2*(h~=0));
Nt = 1024;
[~,~,hinds,rinds,iinds] = HINDS(2,h);

famps = [1e-2 1e-1 1 2 4];

Fl = zeros(2*Nhc, 1);
Fl(rinds(1:2)) = [1;0];
Wst = 0.2;
Wen = 3;
dw = 0.02;
UxwC = cell(size(famps));

epN = 1e0;

% Sopt = struct('stepmax', 500, 'parametrization', 'orthogonal', 'jac', 'full', 'dynamicDscale', 1);
Copt = struct('Nmax', 200, 'DynDscale', 1, 'solverchoice', 1);
for fi=1:length(famps)
    U0 = kron(eye(Nhc), Lst)*(HARMONICSTIFFNESS(Lst'*p.M*Lst, Lst'*p.C*Lst, Lst'*p.K*Lst, Wst, h)\(kron(eye(Nhc),Lst')*Fl*famps(fi)));

    Copt.Dscale = [ones(Nhc*2,1)*famps(fi); 1.0];
    if fi~=1
%         Copt.Dscale = [abs(UxwC{fi-1}(1:end-1,1))*famps(fi)/famps(fi-1)+1e-6; Wst];
        U0 = UxwC{fi-1}(1:end-1,1)*famps(fi)/famps(fi-1);
    end

%     UxwC{fi} = solve_and_continue(U0, @(Uw) HBRESFUN(Uw,Fl*famps(fi),h,Nt,p,epN), Wst, Wen, dw, Sopt);
    UxwC{fi} = CONTINUE(@(Uw) HBRESFUN(Uw,Fl*famps(fi),h,Nt,p,epN), U0, Wst, Wen, dw, Copt);
end
%%

largs = {'-', 'LineWidth', 2};
hi = 1:2:7;

for i=1:length(hi)
    figure((i-1)*10+3)
    clf()
    aa = gobjects(size(famps));
    for fi=1:length(famps)
        try
            subplot(2,1,1)
            plot(UxwC{fi}(end,:), abs(UxwC{fi}(rinds(hi(i)),:)+1j*UxwC{fi}(iinds(hi(i)),:))/famps(fi), largs{:}); hold on
    
            subplot(2,1,2)
            aa(fi) = plot(UxwC{fi}(end,:), rad2deg(angle(UxwC{fi}(rinds(hi(i)),:)+1j*UxwC{fi}(iinds(hi(i)),:))), largs{:}); hold on
            legend(aa(fi), sprintf('F = %.2f N', famps(fi)));
        end
    end
    subplot(2,1,1)
    xlabel('Frequency (rad/s)')
    ylabel('Norm. Amp. (m/N)')
    xlim([Wst Wen])
    grid on
    subplot(2,1,2)
    xlabel('Frequency (rad/s)')
    ylabel('FRF Phase (degs)')
    xlim([Wst Wen])
    grid on
    ll=legend(aa, 'Location', 'best');
    if i~=1
        set(ll, 'Visible', 'off');
    end
    set(gcf, 'Color', 'white')
    if savfig
        export_fig(sprintf('./FIGS/B_FRESP_hi%d_%d.png', hi(i), h(end)), '-dpng');
    end
end
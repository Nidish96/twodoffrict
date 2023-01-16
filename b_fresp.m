clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%% Parameters
bt = 0.1;
nlpars = [0.005;2];  % [alpha, mu]. [1;2]
p = struct('M', eye(2), 'C', bt*[0 0;0 1], 'K', [2 -1;-1 2], ...
    'nlpars', nlpars, 'fv', [1;0]*10, 'epN', 1e2);

[V,Wsr] = eig(p.K, p.M);
[Wsr, si] = sort(sqrt(diag(Wsr)));
V = V(:,si);

%% HB Simulation
h = 0:1;  % Choose harmonics to balance. Ex: h = [0:5];
Nhc = sum((h==0)+2*(h~=0));
Nt = 2048;

famps = [5 10 20];

Fl = zeros(2*Nhc, 1);
Wst = 0.2;
Wen = 4;
dw = 0.1;
UxwC = cell(size(famps));

Sopt = struct('stepmax', 200, 'parametrization', 'orthogonal', 'jac', 'full', 'dynamicDscale', 1);
for fi=1:length(famps)
    Fl(3:4) = [1;0]*famps(fi);

    U0 = HARMONICSTIFFNESS(p.M, p.C, p.K, Wst, h)\Fl;

    p.epN = 1e2;
    UxwC{fi} = solve_and_continue(U0, @(Uw) HBRESFUN(Uw,Fl,h,Nt,p), Wst, Wen, dw, Sopt);
end

figure(1)
clf()
for fi=1:length(famps)
    plot(UxwC{fi}(end,:), abs(UxwC{fi}(3,:)+1j*UxwC{fi}(5,:)), '.-'); hold on
end
%%
rng(1)
Uw = rand(2*Nhc+1,1);
[R, dRdU, dRdw] = HBRESFUN(Uw,Fl,h,Nt,p);

hv = zeros(2*Nhc+1,1);
hm = 1e-6;
Jnum = zeros(2*Nhc, 2*Nhc+1);
for hi=1:2*Nhc+1
    hv(hi) = 1;
    Rp = HBRESFUN(Uw+hv*hm,Fl,h,Nt,p);
    Rm = HBRESFUN(Uw-hv*hm,Fl,h,Nt,p);
    hv(hi) = 0;

    Jnum(:,hi) = (Rp-Rm)/(2*hm);
end
Jan = [dRdU dRdw];
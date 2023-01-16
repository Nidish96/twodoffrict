function [T, U, Ud, Udd, Fnl] = NEWMARKBMARCH(Ft, U0, Ud0, T, p, varargin)
%NEWMARKBMARCH Provides the implicit Newmark-Beta March
%
%   USAGE:
%       [T, U, Ud, Udd] = HHTAMARCH(U0, Ft, T, p, opts);
%   INPUTS:
%       U0  :   (Nd,1)
%       Ft  :   (Nd,Nt)
%       T   :   (1,Nt)
%       p   :
%       opts:
%   OUTPUTS:
%       T   :
%       U   :
%       Ud  :
%       Udd :

% [average acceleration]: beta[1/4], gamma[1/2]
% [implicit linear acceleration]: beta[1/6], gamma[1/2]
% stability condition: 2b >= g >= 1/2

    if isempty(varargin)
        b = 1/4;
        g = 1/2;
    else
        b = varargin{1};
        g = varargin{2};
    end

    Nt = length(T);
    dt = T(2)-T(1);
    Nd = length(U0);
    U  = zeros(Nd, Nt); U(:,1) = U0;
    Ud = zeros(Nd, Nt); Ud(:,1)= Ud0;
    Udd= zeros(Nd, Nt); 
    Fnl= zeros(Nd, Nt);
    % Initial Acceleration
    Fnl(:,1) = [p.nlpars(1)*U0(1)^3; p.nlpars(2)*sign(Ud0(2))];
    Udd(:,1) = p.M\(Ft(:,1)-p.C*Ud(:,1)-p.K*U(:,1)-Fnl(:,1));

    opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
    Uddg = Udd(:,1);
    for ti=2:Nt
        Udd(:,ti) = fsolve(@(udd) TRESFUN(udd, Ft(:,ti), dt, p, b, g, ...
            U(:,ti-1), Ud(:, ti-1), Udd(:, ti-1), Fnl(:, ti-1), Ft(:, ti-1)), ...
        Uddg, opt);

        [~, ~, U(:,ti), Ud(:,ti), Fnl(:,ti)] = TRESFUN(Udd(:,ti), Ft(:,ti), dt, p, b, g, ...
            U(:,ti-1), Ud(:, ti-1), Udd(:, ti-1), Fnl(:, ti-1), Ft(:, ti-1));
        Uddg = Udd(:,ti);
    end
end

% Implicit Res Fun
function [R, dRdUdd, U, Ud, Fnl] = TRESFUN(Udd, Fex, dt, p, b, g, varargin)
    if isempty(varargin)
        Up  = Udd*0;
        Udp = Udd*0;
        Uddp= Udd*0;
        Fnlp= Udd*0;
        Fexp= Udd*0;
    else
        Up  = varargin{1};
        Udp = varargin{2};
        Uddp= varargin{3};
        Fnlp= varargin{4};
        Fexp= varargin{5};
    end

    Z1 = p.M + g*dt*p.C + b*dt^2*p.K;
    Z2 = - p.M + (1-g)*dt*p.C + (0.5-b)*dt^2*p.K;
    Z3 = dt*p.K;
    
    U  = Up + dt*Udp + dt^2/2*((1-2*b)*Uddp+2*b*Udd);
    Ud = Udp+ dt*((1-g)*Uddp+g*Udd);

    Fnl = [p.nlpars(1)*U(1)^3; 0];
    Jnl = [3*p.nlpars(1)*U(1)^2*(dt^2*b) 0;0 0]; 

    R = Z1*Udd + Z2*Uddp + Z3*Udp + (Fnl-Fnlp) - (Fex-Fexp);
    dRdUdd = Z1 + Jnl;

    % Correct friction element
    if abs(-R(2)+p.epN*Ud(2))<p.nlpars(2)
        Fnl(2) = -R(2)+p.epN*Ud(2);
        Jnl(2,:) = -dRdUdd(2,:) + p.epN*[0 1]*dt*g;
    else
        Fnl(2) = p.nlpars(2)*sign(-R(2)+p.epN*Ud(2));
        Jnl(2,:) = 0;
    end

    R = Z1*Udd + Z2*Uddp + Z3*Udp + (Fnl-Fnlp) - (Fex-Fexp);
    dRdUdd = Z1 + Jnl;
end
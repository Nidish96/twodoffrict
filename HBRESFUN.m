function [R, dRdU, dRdw, Fnl] = HBRESFUN(Uw, Fl, h, Nt, p)
%HBRESFUN returns the residue for the forced problem
%
%   USAGE:
%       [R, dRdU, dRdw] = HBRESFUN(Uw, Fl, h, Nt, p);
%   INPUTS:
%       Uw  :
%       Fl  :
%       h   :
%       Nt  :
%       p   :
%   OUTPUTS:
%       R   :
%       dRdU:
%       dRdw:
    
    Nhc = sum((h==0)+2*(h~=0));
    [E, dEdw] = HARMONICSTIFFNESS(p.M, p.C, p.K, Uw(end), h);
    [D1, dD1dw] = HARMONICSTIFFNESS(0, 1, 0, Uw(end), h);

    % Evaluate Nonlinearities
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);

    ut = TIMESERIES_DERIV(Nt, h, reshape(Uw(1:end-1), 2, Nhc)', 0);
    
    % a. Cubic Nonlinearity
    f1    = p.nlpars(1)*ut(:,1).^3;
    dfdu1 = 3*p.nlpars(1)*ut(:,1).^2.*cst;

    F1 = GETFOURIERCOEFF(h, f1);
    dFdU1 = GETFOURIERCOEFF(h, dfdu1);

    % b. Rigid Coulomb Element (predict in f-domain; correct in t-domain)
    I = eye(2*Nhc);
    F2 = Fl(2:2:end)-E(2:2:end,:)*Uw(1:end-1) + p.epN*D1*Uw(2:2:end);
    dFdU2 = -E(2:2:end,:) + p.epN*D1*I(2:2:end,:);
    dFdw2 = -dEdw(2:2:end,:)*Uw(1:end-1) + p.epN*dD1dw*Uw(2:2:end);

    f2 = TIMESERIES_DERIV(Nt, h, F2, 0);
    dfdu2 = TIMESERIES_DERIV(Nt, h, dFdU2, 0);
    dfdw2 = TIMESERIES_DERIV(Nt, h, dFdw2, 0);

    % slip update
    dfdu2(abs(f2)>p.nlpars(2), :) = 0;
    dfdw2(abs(f2)>p.nlpars(2), :) = 0;
    f2(abs(f2)>p.nlpars(2)) = p.nlpars(2)*sign(f2(abs(f2)>p.nlpars(2)));
    
    F2 = GETFOURIERCOEFF(h, f2);
    dFdU2 = GETFOURIERCOEFF(h, dfdu2);
    dFdw2 = GETFOURIERCOEFF(h, dfdw2);

    % Assemble Nonlinear Force
    Fnl = reshape([F1 F2]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, :) = dFdU2;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2;

    % Construct Residual
    R = E*Uw(1:end-1) + Fnl - Fl;
    dRdU = E + Jnl;
    dRdw = dEdw*Uw(1:end-1) + dFnldw;
end
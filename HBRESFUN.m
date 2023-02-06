function [R, dRdU, dRdw, Fnl] = HBRESFUN(Uw, Fl, h, Nt, p, epN)
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
    sct = TIMESERIES_DERIV(Nt, h, D1, 0);

    ut = TIMESERIES_DERIV(Nt, h, reshape(Uw(1:end-1), 2, Nhc)', 0);
    udt = TIMESERIES_DERIV(Nt, h, D1*reshape(Uw(1:end-1), 2, Nhc)', 0);
    duddwt = TIMESERIES_DERIV(Nt, h, dD1dw*reshape(Uw(1:end-1), 2, Nhc)', 0);
    
    % a. Cubic Nonlinearity
    f1    = p.nlpars(1)*ut(:,1).^3;
    dfdu1 = 3*p.nlpars(1)*ut(:,1).^2.*cst;

    F1    = GETFOURIERCOEFF(h, f1);
    dFdU1 = GETFOURIERCOEFF(h, dfdu1);

    % b.1. Cubic damping
    f2    = p.nlpars(2)*udt(:,2).^3;
    dfdu2 = 3*p.nlpars(2)*udt(:,2).^2.*sct;
    dfdw2 = (3*p.nlpars(2)*udt(:,2).^2).*duddwt(:,2);

    F2    = GETFOURIERCOEFF(h, f2);
    dFdU2 = GETFOURIERCOEFF(h, dfdu2);
    dFdw2 = GETFOURIERCOEFF(h, dfdw2);

    % Preliminary Assembly
    Fnl = reshape([F1 F2]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, 2:2:end) = dFdU2;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2;

    % Preliminary Residue
    R = E*Uw(1:end-1)+Fnl-Fl;
    dRdU = E+Jnl;
    dRdw = dEdw*Uw(1:end-1)+dFnldw;

    % b.2. Rigid Coulomb Element (predict in f-domain; correct in t-domain)
    if isinf(p.kt)
        I = eye(2*Nhc);        
        F3 = R(2:2:end) + epN*D1*Uw(2:2:end);
        dFdU3 = dRdU(2:2:end,:) + epN*D1*I(2:2:end,:);
        dFdw3 = dRdw(2:2:end,:) + epN*dD1dw*Uw(2:2:end);
    
        f3 = TIMESERIES_DERIV(Nt, h, F3, 0);
        dfdu3 = TIMESERIES_DERIV(Nt, h, dFdU3, 0);
        dfdw3 = TIMESERIES_DERIV(Nt, h, dFdw3, 0);
    
        % slip update
        dfdu3(abs(f3)>p.nlpars(3), :) = 0;
        dfdw3(abs(f3)>p.nlpars(3), :) = 0;
        f3(abs(f3)>p.nlpars(3)) = p.nlpars(3)*sign(f3(abs(f3)>p.nlpars(3)));
    else
        its = 0;
        f3 = zeros(Nt, 1);
        dfdu3 = zeros(Nt, 2*Nhc);
        dfdw3 = zeros(Nt, 1);
        fprev = 0;
        while its==0 || abs(fprev-f3(1))/p.nlpars(3)<1e-10
            fprev = f3(1);
            for ti=1:Nt
                tim1 = mod(ti-1-1,Nt)+1;

                fsp = p.kt*(ut(ti,2)-ut(tim1,2)) + f3(tim1);
                dfsp = p.kt*(cst(ti,:)-cst(tim1,:)) + dfdu3(tim1,2:2:end);

                if abs(fsp)<p.nlpars(3)
                    f3(ti) = fsp;
                    dfdu3(ti,2:2:end) = dfsp;
                else
                    f3(ti) = p.nlpars(3)*sign(fsp);
                    dfdu3(ti,:) = 0;
                end
            end
            its = its+1;
        end
    end
    F3 = GETFOURIERCOEFF(h, f3);
    dFdU3 = GETFOURIERCOEFF(h, dfdu3);
    dFdw3 = GETFOURIERCOEFF(h, dfdw3);

    % Assemble Nonlinear Force
    Fnl = reshape([F1 F2+F3]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, 2:2:end) = dFdU2;
    Jnl(2:2:end, :) = Jnl(2:2:end, :) + dFdU3;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2+dFdw3;

    % Construct Residual
    R = E*Uw(1:end-1) + Fnl - Fl;
    dRdU = E + Jnl;
    dRdw = dEdw*Uw(1:end-1) + dFnldw;
end
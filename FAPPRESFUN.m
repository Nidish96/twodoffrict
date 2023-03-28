function [R, dRdUw, dRdf, Fnl] = FAPPRESFUN(Uwf, Fl, h, Nt, p, epN)
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
    w = Uwf(end-1);
    lf = Uwf(end);
    f = 10^lf;
    dfdlf = log(10)*f;
    [E, dEdw] = HARMONICSTIFFNESS(p.M, p.C, p.K, w, h);
    [D1, dD1dw] = HARMONICSTIFFNESS(0, 1, 0, w, h);

    U = Uwf(1:end-2)*f;
    dUdf = Uwf(1:end-2)*dfdlf;

    % Evaluate Nonlinearities
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
    sct = TIMESERIES_DERIV(Nt, h, D1, 0);

    ut = TIMESERIES_DERIV(Nt, h, reshape(U, 2, Nhc)', 0);
    udt = TIMESERIES_DERIV(Nt, h, D1*reshape(U, 2, Nhc)', 0);
    duddwt = TIMESERIES_DERIV(Nt, h, dD1dw*reshape(U, 2, Nhc)', 0);
    
    % a. Cubic Nonlinearity
    f1    = p.nlpars(1)*ut(:,1).^3;
    dfdu1 = 3*p.nlpars(1)*ut(:,1).^2.*cst;
    dfdf1 = 3*p.nlpars(1)*ut(:,1).^3/f*dfdlf;

    F1    = GETFOURIERCOEFF(h, f1);
    dFdU1 = GETFOURIERCOEFF(h, dfdu1);
    dFdf1 = GETFOURIERCOEFF(h, dfdf1);

    % b.1. Cubic damping
    f2    = p.nlpars(2)*udt(:,2).^3;
    dfdu2 = 3*p.nlpars(2)*udt(:,2).^2.*sct;
    dfdw2 = (3*p.nlpars(2)*udt(:,2).^2).*duddwt(:,2);
    dfdf2 = 3*p.nlpars(2)*udt(:,2).^3/f*dfdlf;

    F2    = GETFOURIERCOEFF(h, f2);
    dFdU2 = GETFOURIERCOEFF(h, dfdu2);
    dFdw2 = GETFOURIERCOEFF(h, dfdw2);
    dFdf2 = GETFOURIERCOEFF(h, dfdf2);

    % Preliminary Assembly
    Fnl = reshape([F1 F2]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, 2:2:end) = dFdU2;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2;
    dFnldf = zeros(2*Nhc,1);
    dFnldf(2:2:end) = dFdf1+dFdf2;

    % Preliminary Residue
    R = E*U+Fnl-Fl*f;
    dRdU = E+Jnl;
    dRdw = dEdw*U+dFnldw;
    dRdf = E*Uwf(1:end-2)*dfdlf + dFnldf - Fl*dfdlf;

    % b.2. Rigid Coulomb Element (predict in f-domain; correct in t-domain)
    if isinf(p.kt)
        I = eye(2*Nhc);        
        F3 = R(2:2:end) + epN*D1*U(2:2:end);
        dFdU3 = dRdU(2:2:end,:) + epN*D1*I(2:2:end,:);
        dFdw3 = dRdw(2:2:end) + epN*dD1dw*U(2:2:end);
        dFdf3 = dRdf(2:2:end) + epN*D1*dUdf(2:2:end);
    
        f3 = TIMESERIES_DERIV(Nt, h, F3, 0);
        dfdu3 = TIMESERIES_DERIV(Nt, h, dFdU3, 0);
        dfdw3 = TIMESERIES_DERIV(Nt, h, dFdw3, 0);
        dfdf3 = TIMESERIES_DERIV(Nt, h, dFdf3, 0);
    
        % slip update
        dfdu3(abs(f3)>p.nlpars(3), :) = 0;
        dfdw3(abs(f3)>p.nlpars(3), :) = 0;
        dfdf3(abs(f3)>p.nlpars(3), :) = 0;
        f3(abs(f3)>p.nlpars(3)) = p.nlpars(3)*sign(f3(abs(f3)>p.nlpars(3)));
    else
        its = 0;
        f3 = zeros(Nt, 1);
        dfdu3 = zeros(Nt, 2*Nhc);
        dfdw3 = zeros(Nt, 1);
        dfdf3 = zeros(Nt, 1);
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
                    dfdw3(ti) = dfdw3(tim1);
                    dfdf3(ti) = dfdf3(tim1);
                else
                    f3(ti) = p.nlpars(3)*sign(fsp);
                    dfdu3(ti,:) = 0;
                    dfdw3(ti) = 0;
                    dfdf3(ti) = 0;
                end
            end
            its = its+1;
        end
    end
    F3 = GETFOURIERCOEFF(h, f3);
    dFdU3 = GETFOURIERCOEFF(h, dfdu3);
    dFdw3 = GETFOURIERCOEFF(h, dfdw3);
    dFdf3 = GETFOURIERCOEFF(h, dfdf3);

    % Assemble Nonlinear Force
    Fnl = reshape([F1 F2+F3]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, 2:2:end) = dFdU2;
    Jnl(2:2:end, :) = Jnl(2:2:end, :) + dFdU3;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2+dFdw3;
    dFnldf = zeros(2*Nhc,1);
    dFnldf(1:2:end) = dFdf1;
    dFnldf(2:2:end) = dFdf2+dFdf3;

    % Construct Residual
    R = [E*U+Fnl-Fl*f;
        Fl'*U];
    dRdUw = [(E+Jnl)*f dEdw*U+dFnldw;
        Fl'*f 0];
    dRdf = [E*Uwf(1:end-2)*dfdlf+dFnldf-Fl*dfdlf;
        Fl'*dUdf];
end
function [R, dRdUxw, dRdlq, Fnl] = EPMCRESFUN(Uxwq, Fl, h, Nt, p, epN)
%HBRESFUN returns the residue for the forced problem
%
%   USAGE:
%       [R, dRdU, dRdw] = EPMCRESFUN(Uxwq, Fl, h, Nt, p);
%   INPUTS:
%       Uxwq:
%       Fl  :
%       h   :
%       Nt  :
%       p   :
%   OUTPUTS:
%       R   :
%       dRdU:
%       dRdw:
    
    Nhc = sum((h==0)+2*(h~=0));
    
    lq = Uxwq(end);
    q = 10^lq;
    dqdlq = log(10)*q;
    
    [~, ~, hinds, rinds, iinds] = HINDS(2, h(:));
    Asm = zeros(Nhc*2,1);  Asc = zeros(Nhc*2,1);
    Asm(hinds) = 1;
    Asc([rinds iinds]) = 1;

    Usc = Uxwq(1:end-3).*(Asm+Asc*q);

    w = Uxwq(end-1);
    xi = Uxwq(end-2);

    [E, dEdw] = HARMONICSTIFFNESS(p.M, p.C-xi*p.M, p.K, w, h);
    dEdxi = HARMONICSTIFFNESS(zeros(2), -p.M, zeros(2), w, h);
    [D1, dD1dw] = HARMONICSTIFFNESS(0, 1, 0, w, h);

    % Evaluate Nonlinearities
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
    sct = TIMESERIES_DERIV(Nt, h, D1, 0);

    ut = TIMESERIES_DERIV(Nt, h, reshape(Usc, 2, Nhc)', 0);
    udt = TIMESERIES_DERIV(Nt, h, D1*reshape(Usc, 2, Nhc)', 0);
    duddwt = TIMESERIES_DERIV(Nt, h, dD1dw*reshape(Usc, 2, Nhc)', 0);
    
    % a. Cubic Nonlinearity
    f1    = p.nlpars(1)*ut(:,1).^3;
    dfdu1 = 3*p.nlpars(1)*ut(:,1).^2.*cst;

    F1 = GETFOURIERCOEFF(h, f1);
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
    R = E*Usc+Fnl;
    dRdUxw = [(E+Jnl)*diag(Asm+Asc*q) dEdxi*Usc dEdw*Usc+dFnldw];
    dRdlq  = (E+Jnl)*(Uxwq(1:end-3).*Asc*dqdlq);

    if isinf(p.kt)
        % b. Rigid Coulomb Element (predict in f-domain; correct in t-domain)
        I = eye(2*Nhc);
        F3 = -R(2:2:end) + epN*D1*Usc(2:2:end);
        dFdUxw3 = -dRdUxw(2:2:end, :) + epN*[D1*I(2:2:end,:)*diag(Asm+Asc*q) zeros(Nhc,1) dD1dw*Usc(2:2:end)];
        dFdlq = -dRdlq(2:2:end) + epN*D1*(Uxwq(2:2:end-3).*Asc(2:2:end)*dqdlq);

        f3 = TIMESERIES_DERIV(Nt, h, F3, 0);
        dfdUxw3 = TIMESERIES_DERIV(Nt, h, dFdUxw3, 0);
        dfdlq3 = TIMESERIES_DERIV(Nt, h, dFdlq, 0);
    
        % Slip update
        dfdUxw3(abs(f3)>=p.nlpars(3), :) = 0;
        dfdlq3(abs(f3)>=p.nlpars(3), :) = 0;
        f3(abs(f3)>=p.nlpars(3)) = p.nlpars(3)*sign(f3(abs(f3)>=p.nlpars(3)));

        F3 = GETFOURIERCOEFF(h, f3);
        dFdUxw3 = GETFOURIERCOEFF(h, dfdUxw3);
        dFdlq3 = GETFOURIERCOEFF(h, dfdlq3);
    else  % Elastic Dry Friction Element
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
        dfdxi3 = zeros(Nt,1);

        F3 = GETFOURIERCOEFF(h, f3);
        dFdUxw3 = GETFOURIERCOEFF(h, [dfdu3*diag(Asm+Asc*q) dfdxi3 dfdw3]);
        dFdlq3 = GETFOURIERCOEFF(h, dfdu3)*(Uxwq(2:2:end-3).*Asc(2:2:end)*dqdlq);    
    end

    % Assemble Nonlinear Force
    Fnl = reshape([F1 F2+F3]', 2*Nhc, 1);
    dFnldUxw = zeros(2*Nhc, 2*Nhc+2);
    dFnldUxw(1:2:end, 1:2:2*Nhc) = dFdU1;
    dFnldUxw(2:2:end, [2:2:2*Nhc end]) = [dFdU2 dFdw2];
    dFnldUxw(2:2:end, :) = dFnldUxw(2:2:end, :) + dFdUxw3;
    dFnldlq = zeros(2*Nhc,1);
    dFnldlq(1:2:end) = dFdU1*(Uxwq(1:2:end-3).*Asc(1:2:end)*dqdlq);
    dFnldlq(2:2:end) = dFdU2*(Uxwq(2:2:end-3).*Asc(2:2:end)*dqdlq) + dFdlq3;

    % First harmonics Amplitude Constraint
    acons = Uxwq(rinds(1:2))'*p.M*Uxwq(rinds(1:2)) + Uxwq(iinds(1:2))'*p.M*Uxwq(iinds(1:2)) - 1;
    daconsdU = zeros(1, Nhc*2+2);
    daconsdU(rinds(1:2)) = 2*Uxwq(rinds(1:2))'*p.M;
    daconsdU(iinds(1:2)) = 2*Uxwq(iinds(1:2))'*p.M;
    daconsdlq = 0;

    % Construct Residual
    R = [E*Usc+Fnl;
        acons;
        Fl'*Uxwq(1:end-3)];
    dRdUxw = [[E*diag(Asm+Asc*q) dEdxi*Usc dEdw*Usc]+dFnldUxw;
        daconsdU;
        Fl' 0 0];
    dRdlq = [E*(Uxwq(1:end-3).*Asc*dqdlq)+dFnldlq;
        daconsdlq;
        0];
end
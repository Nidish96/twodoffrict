function [R, dRdUxw, dRdlq, Fnl] = EPMCRESFUN(Uxwq, Fl, h, Nt, p)
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
    
    hi = find(h==1);
    H1sel = eye(2*Nhc);
    if h(1)==0
        Asm = [ones(2,1); zeros(2*(Nhc-1),1)];
        Asc = [zeros(2,1); ones(2*(Nhc-1),1)];

        H1sel = H1sel(2 + (hi-2)*4+(1:4), :);
    else
        Asm = zeros(2*Nhc,1);
        Asc = ones(2*Nhc,1);

        H1sel = H1sel((hi-1)*4+(1:4), :);
    end
    Usc = Uxwq(1:end-3).*(Asm+Asc*q);

    w = Uxwq(end-1);
    xi = Uxwq(end-2);

    [E, dEdw] = HARMONICSTIFFNESS(p.M, p.C-xi*p.M, p.K, w, h);
    dEdxi = HARMONICSTIFFNESS(p.M*0, -p.M, p.K*0, w, h);
    [D1, dD1dw] = HARMONICSTIFFNESS(0, 1, 0, w, h);

    % Evaluate Nonlinearities
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);

    ut = TIMESERIES_DERIV(Nt, h, reshape(Usc, 2, Nhc)', 0);
    
    % a. Cubic Nonlinearity
    f1    = p.nlpars(1)*ut(:,1).^3;
    dfdu1 = 3*p.nlpars(1)*ut(:,1).^2.*cst;

    F1 = GETFOURIERCOEFF(h, f1);
    dFdU1 = GETFOURIERCOEFF(h, dfdu1);

    % b. Rigid Coulomb Element (predict in f-domain; correct in t-domain)
    I = eye(2*Nhc);
    F2 = Fl(2:2:end)-E(2:2:end,:)*Usc + p.epN*D1*Usc(2:2:end);
    dFdU2 = -E(2:2:end,:) + p.epN*D1*I(2:2:end,:);
    dFdw2 = -dEdw(2:2:end,:)*Usc + p.epN*dD1dw*Usc(2:2:end);
    dFdxi2 = -dEdxi(2:2:end,:)*Usc;

    f2 = TIMESERIES_DERIV(Nt, h, F2, 0);
    dfdu2 = TIMESERIES_DERIV(Nt, h, dFdU2, 0);
    dfdw2 = TIMESERIES_DERIV(Nt, h, dFdw2, 0);
    dfdxi2 = TIMESERIES_DERIV(Nt, h, dFdxi2, 0);

    % slip update
    dfdu2(abs(f2)>p.nlpars(2), :) = 0;
    dfdw2(abs(f2)>p.nlpars(2), :) = 0;
    dfdxi2(abs(f2)>p.nlpars(2), :) = 0;
    f2(abs(f2)>p.nlpars(2)) = p.nlpars(2)*sign(f2(abs(f2)>p.nlpars(2)));
    
    F2 = GETFOURIERCOEFF(h, f2);
    dFdU2 = GETFOURIERCOEFF(h, dfdu2);
    dFdw2 = GETFOURIERCOEFF(h, dfdw2);
    dFdxi2 = GETFOURIERCOEFF(h, dfdxi2);

    % Assemble Nonlinear Force
    Fnl = reshape([F1 F2]', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = dFdU1;
    Jnl(2:2:end, :) = dFdU2;
    dFnldw = zeros(2*Nhc,1);
    dFnldw(2:2:end) = dFdw2;
    dFnldxi= zeros(2*Nhc,1);
    dFnldxi(2:2:end) = dFdxi2;
    
    % Construct Residual
    R = [E*Usc+Fnl;
        (H1sel*Usc)'*blkdiag(p.M,p.M)*(H1sel*Usc)-q^2;
        Fl'*Uxwq(1:end-3)];
    dRdUxw = [(E+Jnl).*(Asm+Asc*q) dEdxi*Usc+dFnldxi dEdw*Usc+dFnldw;
        2*(H1sel*Usc)'*blkdiag(p.M,p.M)*(H1sel.*(Asm+Asc*q)') 0 0;
        Fl' 0 0];
    dRdlq = [(E+Jnl)*(Uxwq(1:end-3).*Asc*dqdlq);
        2*(H1sel*Usc)'*blkdiag(p.M,p.M)*(H1sel*(Uxwq(1:end-3).*Asc*dqdlq))-2*q*dqdlq;
        0];
end
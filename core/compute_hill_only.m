function Hill = compute_hill_only(p, CeP, CeR, model, ~)
    epsv = 1e-12;
    C50P = p(1); C50R = p(2); gamma = p(3); s4 = p(4);
    CeP_mgL = max(CeP, 0);
    CeR_ngmL = 1000 * max(CeR, 0);
    
    if strcmpi(model, 'vanluchene')
        beta = s4;
        Ip = CeP_mgL / max(C50P, epsv);
        Ir = CeR_ngmL / max(C50R, epsv);
        denom = Ip + Ir + epsv;
        theta = max(0, min(1, Ip / denom));
        w = 1 - beta * theta + beta * (theta.^2) + epsv;
        U = (Ip + Ir) ./ w;
    else
        alpha = s4;
        termP = CeP_mgL / max(C50P, epsv);
        termR = CeR_ngmL / max(C50R, epsv);
        U = termP + termR + alpha * termP * termR + epsv;
    end
    
    Hill = U.^gamma ./ (1 + U.^gamma);
end

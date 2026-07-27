function [bis, H, xP, xR] = predict_bis_loglin3d(p, CeP, CeR, E0, BISmin, state)
% H, xP and xR are returned so the EKF builds its Jacobian from this evaluation.

    a0 = p(1);
    aP = p(2);
    aR = p(3);

    xP = log(CeP + state.eps_P);
    xR = log(CeR * 1000 + state.eps_R);

    Z = a0 + aP * xP + aR * xR;
    Z = max(-8, min(8, Z));
    H = 1 / (1 + exp(-Z));

    bis = BISmin + (E0 - BISmin) * (1 - H);
    bis = max(0, min(100, bis));
end

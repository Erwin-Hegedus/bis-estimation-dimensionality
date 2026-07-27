function S = param_scales()
%PARAM_SCALES  Per-model parameter scales used to build the process noise.
%   S = PARAM_SCALES() returns the covariance scales S.m1d .. S.m4d, so that a
%   sweep over a single scalar q sets each model's process noise to q * S. The
%   scales are the initial parameter standard deviations of each estimator
%   (init_ekf_kscale, init_ekf_2d, init_ekf_loglin3d, init_ekf_4d), which
%   makes q dimensionless and comparable across models of different dimension.
%
%   S.sd_* are the same quantities as vectors, for normalizing parameter drift.

    S.sd_1d = 0.20;
    S.sd_2d = [0.3; 0.3];
    S.sd_3d = [1.0; 0.5; 0.1];
    S.sd_4d = [1.0; 1.5; 0.4; 0.15];

    S.m1d = S.sd_1d^2;
    S.m2d = diag(S.sd_2d.^2);
    S.m3d = diag(S.sd_3d.^2);
    S.m4d = diag(S.sd_4d.^2);
end

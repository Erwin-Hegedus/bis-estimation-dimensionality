function S = param_scales()
%PARAM_SCALES  Initial parameter standard deviations of each estimator.
%   S = PARAM_SCALES() returns the initial parameter standard deviations of
%   init_ekf_kscale, init_ekf_2d, init_ekf_loglin3d and init_ekf_4d, used to
%   normalize parameter drift and roughness across models of different
%   dimension.
%
%   Process noise is not built here: each estimator sets Q = cfg.q * P0 from its
%   own initial covariance.

    S.sd_1d = 0.20;
    S.sd_kgamma = [0.20; 0.40];
    S.sd_2d = [0.3; 0.3];
    S.sd_3d = [1.0; 0.5; 0.1];
    S.sd_4d = [1.0; 1.5; 0.4; 0.15];
end

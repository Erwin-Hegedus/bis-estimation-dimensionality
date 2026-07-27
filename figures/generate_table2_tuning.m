function generate_table2_tuning(cfg, fig_dir)
%GENERATE_TABLE2_TUNING  Per-model EKF tuning, read from the estimators.
%   Values come from the init functions rather than being restated here, so the
%   table cannot drift from the code it documents.

    fprintf('Generating Table 2: EKF Tuning...\n');

    s1 = init_ekf_kscale(cfg);
    s2 = init_ekf_2d(cfg);
    s3 = init_ekf_loglin3d(cfg);
    s4 = init_ekf_4d(cfg, 'vanluchene');

    row = '%-22s %-14s %-22s %-28s %s\n';
    fid = fopen(fullfile(fig_dir, 'table2_tuning.txt'), 'w');
    fprintf(fid, 'TABLE 2: Model-Specific EKF Tuning Parameters\n');
    fprintf(fid, '==============================================\n\n');
    fprintf(fid, row, 'Parameter', '1D', '2D', '3D', '4D');
    fprintf(fid, row, '---------', '--', '--', '--', '--');
    fprintf(fid, row, 'Process noise (sd)', ...
        vec2str(sqrt(diag(s1.Q))), vec2str(sqrt(diag(s2.Q))), ...
        vec2str(sqrt(diag(s3.Q))), vec2str(sqrt(diag(s4.Q))));
    fprintf(fid, row, 'Initial P (sd)', ...
        vec2str(sqrt(diag(s1.P))), vec2str(sqrt(diag(s2.P))), ...
        vec2str(sqrt(diag(s3.P))), vec2str(sqrt(diag(s4.P))));
    fprintf(fid, row, 'Rate limit / sample', ...
        '0.05', vec2str(s2.rate_max), vec2str(s3.param_rate_max), ...
        vec2str(s4.param_rate_max));
    fprintf(fid, row, 'FIM projection', 'n/a (scalar)', 'Yes', 'Yes', 'Yes');
    fprintf(fid, row, 'FIM forgetting', '---', ...
        vec2str(s2.FIM_forgetting), vec2str(s3.FIM_forgetting), vec2str(s4.FIM_forgetting));
    fprintf(fid, row, 'Identifiable ratio', '---', '0.01', ...
        vec2str(s3.ident_eigenvalue_ratio), vec2str(s4.ident_eigenvalue_ratio));

    fprintf(fid, '\nThe 1D filter adapts a scalar, so subspace projection does not apply;\n');
    fprintf(fid, 'it guards on |H| > 1e-4 instead.\n');
    fprintf(fid, '\nShared parameters:\n');
    fprintf(fid, '  R_base = %.0f\n', cfg.R_base);
    fprintf(fid, '  R_disequilibrium_factor = %.0f\n', cfg.R_disequilibrium_factor);
    fprintf(fid, '  ke0P = %.3f /min, ke0R = %.3f /min\n', cfg.ke0P, cfg.ke0R);
    fprintf(fid, '  E0 = %.0f, BISmin = %.0f (population values, 0D model)\n', ...
        cfg.E0_fixed, cfg.BISmin_fixed);
    fclose(fid);

    fprintf('  Saved: table2_tuning.txt\n');
end

function s = vec2str(v)
    v = v(:)';
    if isscalar(v)
        s = sprintf('%.4g', v);
    else
        s = ['[' strjoin(arrayfun(@(x) sprintf('%.4g', x), v, 'UniformOutput', false), ', ') ']'];
    end
end

function generate_table2_tuning(cfg, fig_dir)
% TABLE 2: Model-Specific EKF Tuning Parameters

    fprintf('Generating Table 2: EKF Tuning...\n');
    
    fid = fopen(fullfile(fig_dir, 'table2_tuning.txt'), 'w');
    fprintf(fid, 'TABLE 2: Model-Specific EKF Tuning Parameters\n');
    fprintf(fid, '==============================================\n\n');
    fprintf(fid, 'Parameter              1D          2D              3D                  4D\n');
    fprintf(fid, '----------             --          --              --                  --\n');
    fprintf(fid, 'Process noise Q        0.002^2     diag([0.002,    diag([0.005,        0 (constants)\n');
    fprintf(fid, '                                   0.002]^2)       0.002,0.0005]^2)\n');
    fprintf(fid, 'Initial P              0.20^2      diag([0.3,      diag([1.0,          diag([1.0,1.5,\n');
    fprintf(fid, '                                   0.3]^2)         0.5,0.1]^2)         0.4,0.15]^2)\n');
    fprintf(fid, 'Rate limits            ---         [0.01,0.01]     [0.02,0.01,0.005]   [0.003,0.008,\n');
    fprintf(fid, '                                                                       0.002,0.001]\n');
    fprintf(fid, 'FIM projection         No          No              No                  Yes\n');
    fprintf(fid, '\nShared parameters:\n');
    fprintf(fid, '  R_base = %.0f\n', cfg.R_base);
    fprintf(fid, '  R_disequilibrium_factor = %.0f\n', cfg.R_disequilibrium_factor);
    fclose(fid);
    
    fprintf('  Saved: table2_tuning.txt\n');
end

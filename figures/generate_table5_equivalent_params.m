function generate_table5_equivalent_params(eq_info, pid, fig_dir, n_examples)
% TABLE 5: Example equivalent 4D parameter sets (tab:equivalent_params)

    if nargin < 4, n_examples = 5; end

    fprintf('Generating Table 5: Equivalent Parameter Sets...\n');

    fid = fopen(fullfile(fig_dir, 'table5_equivalent_params.txt'), 'w');
    fprintf(fid, 'TABLE 5: Example Equivalent Parameter Sets (Patient %d)\n', pid);
    fprintf(fid, '============================================================\n\n');

    if ~eq_info.has_eq || isempty(eq_info.theta_good)
        fprintf(fid, 'No equivalent parameter sets found within %.0f%% of MAE_ref.\n', 10);
        fclose(fid);
        fprintf('  Saved: table5_equivalent_params.txt (empty)\n');
        return;
    end

    n = min(n_examples, eq_info.N_eq);

    fprintf(fid, 'Reference MAE      %.2f BIS\n', eq_info.mae_ref);
    fprintf(fid, 'Equivalent sets    %d / %d (%.1f%%)\n', ...
        eq_info.N_eq, eq_info.N_total, 100 * eq_info.N_eq / eq_info.N_total);
    fprintf(fid, 'Radius (L2)        mean %.2f   max %.2f\n\n', ...
        eq_info.radius_mean, eq_info.radius_max);

    fprintf(fid, '%-14s %-14s %-11s %-11s\n', 'C50P (ug/mL)', 'C50R (ng/mL)', 'gamma', 'beta');
    fprintf(fid, '%-14s %-14s %-11s %-11s\n', '------------', '------------', '-----', '----');
    ref = arrayfun(@(v) sprintf('%.2f (ref)', v), eq_info.theta_ref, 'UniformOutput', false);
    fprintf(fid, '%-14s %-14s %-11s %-11s\n', ref{:});
    for i = 1:n
        row = arrayfun(@(v) sprintf('%.2f', v), eq_info.theta_good(i, :), 'UniformOutput', false);
        fprintf(fid, '%-14s %-14s %-11s %-11s\n', row{:});
    end

    fclose(fid);
    fprintf('  Saved: table5_equivalent_params.txt\n');
end

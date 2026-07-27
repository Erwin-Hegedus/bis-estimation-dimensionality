%RUN_Q_SWEEP_4D  4D Bouillon accuracy and parameter drift against process noise.
%   Sweeps the process-noise scale q and reports cohort MAE and drift, to test
%   whether the 4D result depends on it having been run frozen (Q = 0).

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
cfg        = P.cfg;
DATA_FILE  = P.data_file;
valid_pids = P.cohort;
review_dir = P.outdir;

Q_scale = getfield(param_scales(), 'm4d');
q_grid  = [0, 1e-3, 1e-2, 1e-1, 1, 10];

N = get_patient_count(DATA_FILE);
nq = numel(q_grid);
mae   = nan(N, nq);
drift = nan(N, nq);

for j = 1:nq
    Qj = q_grid(j) * Q_scale;
    fprintf('\n===== Q sweep point %d/%d: q = %.0e =====\n', j, nq, q_grid(j));
    t0 = tic;
    [mae(:,j), drift(:,j)] = run_cohort_van(cfg, Qj, DATA_FILE, N);
    fprintf('q = %.0e : cohort 4D MAE mean %.3f median %.3f (N=%d) in %.0f s\n', ...
        q_grid(j), mean(mae(valid_pids,j),'omitnan'), median(mae(valid_pids,j),'omitnan'), ...
        sum(~isnan(mae(valid_pids,j))), toc(t0));

    sweep = struct('q_grid', q_grid, 'Q_scale', Q_scale, ...
                   'valid_pids', valid_pids, 'mae', mae, 'drift', drift);
    prov = provenance_stamp(cfg, DATA_FILE);
    save(fullfile(review_dir, 'q_sweep_4d.mat'), 'sweep', 'prov');
end

mae_med  = median(mae(valid_pids,:), 1, 'omitnan');
mae_mean = mean(mae(valid_pids,:), 1, 'omitnan');
fprintf('\n=== 4D MAE vs q (cohort N=209) ===\n');
for j = 1:nq
    fprintf('  q = %.0e :  mean %.3f   median %.3f\n', q_grid(j), mae_mean(j), mae_med(j));
end

fig = figure('Color','w','Position',[80 80 1100 460]);

subplot(1,2,1);
semilogx(max(q_grid,1e-4), mae_mean, '-o', 'LineWidth', 2, 'Color', [0.2 0.4 0.7]);
hold on;
xlabel('Process noise scale q', 'FontSize', 14);
ylabel('4D cohort MAE (BIS)', 'FontSize', 14);
title('(a) 4D accuracy vs adaptation', 'FontSize', 15, 'FontWeight', 'bold');
legend({'4D MAE'}, 'Location','best','FontSize',11);
set(gca,'FontSize',12); grid on;

subplot(1,2,2);
semilogx(max(q_grid,1e-4), median(drift(valid_pids,:),1,'omitnan'), '-o', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
xlabel('Process noise scale q', 'FontSize', 14);
ylabel('4D parameter drift (normalized)', 'FontSize', 14);
title('(b) Parameter drift vs adaptation', 'FontSize', 15, 'FontWeight', 'bold');
set(gca,'FontSize',12); grid on;

sgtitle('4D Q-sweep (N=209)', 'FontSize', 16, 'FontWeight', 'bold');
export_figure_png(fig, fullfile(review_dir, 'figure_q_sweep_4d.png'));
fprintf('Saved: %s and figure_q_sweep_4d.png\n', fullfile(review_dir, 'q_sweep_4d.mat'));
